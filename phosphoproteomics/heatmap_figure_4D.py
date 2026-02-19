"""
Figure 4D -- Per-module phosphoproteomic heatmaps (raw phospho).

Replicates the R script phospho_heatmaps_figure_4D.qmd in Python.

Generates per-module heatmaps of curated Figure 4D gene set:
  - 88 curated genes assigned to 8 primary modules
  - Raw phospho (sample == "phospho"), per-timepoint DMSO-normalized log2FC
    (log2(Tryptamine_Xmin / DMSO_Xmin))
  - Median across replicates per timepoint
  - Duplicate phosphopeptides at same site merged by median
  - Rows ordered by log2FC at 15 min (descending)

Required files (relative to repo root):
  1. phosphoproteomics/data/limma_results_annotated.tsv
  2. phosphoproteomics/data/mdata.rds  (read via pyreadr)
  3. phosphoproteomics/data/figure_4D_primary_module_assignments.tsv
  4. data/spongilla_gene_names_final.tsv

Usage:
  conda env create -f environment.yml && conda activate monoamine-sponges
  python phosphoproteomics/heatmap_figure_4D.py
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyreadr
import seaborn as sns

# ============================================================================
# Configuration
# ============================================================================

SAMPLE_TYPE = "phospho"  # raw phospho
MAX_ABS_CAP = 2.1        # cap symmetric color scale at +/- 2.1
CELL_SIZE_INCHES = 0.095  # each square is 0.095 inch x 0.095 inch
FONT_SIZE = 7            # title font size
FONT_SIZE_LABELS = 6     # protein name labels font size

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE,
    'ytick.labelsize': FONT_SIZE,
    'svg.fonttype': 'none',  # keep text as text in SVG (not paths)
})

# The 5 key comparisons for hit/candidate identification (comparison column)
RELEVANT_COMPARISONS = [
    "Tryptamine_3min - DMSO_3min",
    "Tryptamine_15min - DMSO_15min",
    "Tryptamine_30min - DMSO_30min",
    "(Tryptamine_15min - DMSO_15min) - (Tryptamine_3min - DMSO_3min)",
    "(Tryptamine_30min - DMSO_30min) - (Tryptamine_15min - DMSO_15min)",
]

DPI = 300

# ============================================================================
# Paths
# ============================================================================

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent

dir_data = _SCRIPT_DIR / "data"
dir_annot = _REPO_ROOT / "data"
dir_out = _REPO_ROOT / "outs" / "phosphoproteomics" / "fig4d_heatmaps"

dir_out.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Utility functions
# ============================================================================

def normalize_name(x):
    """Normalize automated_name / seurat_name for join matching.

    Strips apostrophes (unicode U+2019 and ASCII) and converts underscore
    format (c12345_g1) to dash format (c12345-g1) for compatibility between
    older automated_name and current seurat_name columns.
    """
    if pd.isna(x):
        return x
    x = str(x).replace("\u2019", "").replace("'", "")
    # Convert underscore format to dash format in Trinity ID prefix
    x = re.sub(r"^(c\d+)_(g\d+)", r"\1-\2", x)
    return x


def slugify(text):
    """Convert text to filename-safe slug."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


# ============================================================================
# Data loading
# ============================================================================

def load_data():
    """Load all required datasets and return processed DataFrames."""

    # 1. Figure 4D module assignments (curated)
    path_assignments = dir_data / "figure_4D_primary_module_assignments.tsv"

    assignments = pd.read_csv(path_assignments, sep="\t")
    assignments = assignments[["automated_name", "Gene.short", "suggested_primary"]].copy()
    assignments = assignments.rename(columns={"suggested_primary": "primary_module"})
    assignments["automated_name_norm"] = assignments["automated_name"].apply(normalize_name)

    print(f"Figure 4D gene assignments: {len(assignments)} genes")
    print(f"Primary modules: {assignments['primary_module'].nunique()}")
    print("\nGenes per module:")
    print(assignments.groupby("primary_module").size().sort_values(ascending=False).to_string())
    print()

    # 2. Limma results
    path_res = dir_data / "limma_results_annotated.tsv"

    print("Loading limma results...")
    res = pd.read_csv(path_res, sep="\t", low_memory=False)

    # 3. Gene annotations (final gene names with paralog suffixes)
    # Join on automated_name (= seurat_name) for robust gene-level matching
    path_gene_annotations = dir_annot / "spongilla_gene_names_final.tsv"

    gene_annotations = pd.read_csv(
        path_gene_annotations, sep="\t",
        usecols=["seurat_name", "Zang_et_al_2026", "name_type"],
    ).rename(columns={
        "Zang_et_al_2026": "Gene_short_new",
        "name_type": "name_type_new",
    })
    gene_annotations["automated_name_norm"] = gene_annotations["seurat_name"].apply(normalize_name)
    gene_annotations = gene_annotations.drop_duplicates(subset="automated_name_norm")
    gene_annotations = gene_annotations[["automated_name_norm", "Gene_short_new", "name_type_new"]]

    # Update Gene.short from final gene names
    res["automated_name_norm"] = res["automated_name"].apply(normalize_name)
    res = res.merge(gene_annotations, on="automated_name_norm", how="left")
    res["Gene.short"] = res["Gene_short_new"].fillna(res["Gene.short"])
    if "name_type" in res.columns:
        res["name_type"] = res["name_type_new"].fillna(res["name_type"])
    else:
        res["name_type"] = res["name_type_new"].fillna("trinity_only")
    res = res.drop(columns=["Gene_short_new", "name_type_new"])

    # 4. Measurement data (mdata) — read from RDS (compact) via pyreadr
    path_mdata = dir_data / "mdata.rds"

    print("Loading mdata from RDS (this may take a moment)...")
    mdata = pyreadr.read_r(str(path_mdata))[None]

    print(f"\nData loaded:")
    print(f"  res: {len(res)} rows")
    print(f"  mdata: {len(mdata)} rows")

    return assignments, res, mdata


# ============================================================================
# Processing pipeline
# ============================================================================

def identify_significant_sites(res):
    """Find all significant phosphosites across the 5 key comparisons."""
    mask = (
        (res["sample"] == SAMPLE_TYPE)
        & (res["comparison"].isin(RELEVANT_COMPARISONS))
        & (res["hit_annotation"].isin(["hit", "candidate"]))
    )
    significant_sites = res.loc[mask, "sequence.id"].unique()
    print(f"\nAll significant phosphosites (raw phospho): {len(significant_sites)}")
    return set(significant_sites)


def build_gene_mapping(res):
    """Build sequence.id -> gene annotation mapping for raw phospho."""
    mask = res["sample"] == SAMPLE_TYPE
    cols = ["sequence.id", "Gene", "Gene.short", "automated_name", "name_type",
            "phospho.position"]
    mapping = res.loc[mask, cols].copy()
    mapping["automated_name_norm"] = mapping["automated_name"].apply(normalize_name)
    mapping = mapping.drop_duplicates()
    return mapping


def filter_figure_4d_sites(gene_mapping, significant_sites, assignments):
    """Filter to Figure 4D genes AND significant sites."""
    # Keep only significant sites
    in_sig = gene_mapping["sequence.id"].isin(significant_sites)
    gene_mapping_sig = gene_mapping[in_sig].copy()

    # Inner join with assignments on normalized automated_name
    sites_fig4d = gene_mapping_sig.merge(
        assignments[["automated_name_norm", "primary_module"]],
        on="automated_name_norm",
        how="inner",
    )

    print(f"Figure 4D phosphosites (hit/candidate): {sites_fig4d['sequence.id'].nunique()}")
    print(f"Figure 4D unique genes: {sites_fig4d['Gene.short'].nunique()}")
    print("\nPhosphosites per module:")
    print(sites_fig4d.groupby("primary_module")["sequence.id"].nunique()
          .sort_values(ascending=False).to_string())
    print()

    return sites_fig4d


def compute_log2fc_matrix(mdata, sites_fig4d):
    """Compute per-timepoint DMSO-normalized log2FC from mdata.

    Each timepoint is normalized to its matched DMSO control:
    log2FC = log2(Tryptamine_Xmin / DMSO_Xmin), computed by subtracting
    median log2(ctrl.ratio) of DMSO replicates from Tryptamine replicates
    at each timepoint (the shared DMSO_3min baseline cancels out).
    """
    site_ids = set(sites_fig4d["sequence.id"].unique())

    # Filter mdata: raw phospho, ctrl.ratio, ALL conditions (Tryptamine + DMSO)
    mask = (
        (mdata["sample"] == SAMPLE_TYPE)
        & (mdata["measurement"] == "ctrl.ratio")
        & (mdata["sequence.id"].isin(site_ids))
    )
    tc_data = mdata[mask].copy()

    # Extract timepoint and treatment
    tc_data["time_min"] = tc_data["condition"].str.extract(r"(\d+)").astype(int)
    tc_data["treatment"] = np.where(
        tc_data["condition"].str.startswith("Tryptamine"), "Tryptamine", "DMSO"
    )
    tc_data["log2_ctrl_ratio"] = np.log2(tc_data["value"])

    # Remove NA and infinite values
    tc_data = tc_data[tc_data["log2_ctrl_ratio"].notna() & np.isfinite(tc_data["log2_ctrl_ratio"])]

    print(f"Timecourse data: {len(tc_data)} observations")

    # Step 1: Median across replicates per sequence.id per treatment per timepoint
    tc_median = (
        tc_data
        .groupby(["sequence.id", "time_min", "treatment"])["log2_ctrl_ratio"]
        .median()
        .reset_index()
    )

    # Step 2: Pivot treatment to columns for per-timepoint subtraction
    tc_wide_treatment = tc_median.pivot_table(
        index=["sequence.id", "time_min"],
        columns="treatment",
        values="log2_ctrl_ratio",
    ).reset_index()
    tc_wide_treatment.columns.name = None

    # Step 3: Per-timepoint DMSO normalization
    # log2(Tryp_Xmin / DMSO_Xmin) = log2(Tryp/DMSO_3min) - log2(DMSO_matched/DMSO_3min)
    tc_wide_treatment["log2FC"] = tc_wide_treatment["Tryptamine"] - tc_wide_treatment["DMSO"]

    # Step 4: Pivot to wide format
    log2fc_wide = tc_wide_treatment[["sequence.id", "time_min", "log2FC"]].pivot(
        index="sequence.id",
        columns="time_min",
        values="log2FC",
    ).reset_index()

    log2fc_wide.columns.name = None
    log2fc_wide = log2fc_wide.rename(columns={
        3: "log2FC_3",
        15: "log2FC_15",
        30: "log2FC_30",
    })

    print(f"Before merging peptides: {len(log2fc_wide)} peptides x 3 timepoints")

    return log2fc_wide


def merge_and_label(log2fc_wide, sites_fig4d):
    """Merge duplicate phosphopeptides at same site and create row labels."""

    # Join with gene/module annotations (still at peptide level)
    join_cols = ["sequence.id", "Gene.short", "primary_module", "name_type", "phospho.position"]
    # Deduplicate sites_fig4d to avoid many-to-many
    sites_dedup = sites_fig4d[join_cols].drop_duplicates(subset="sequence.id")

    heatmap_peptide = log2fc_wide.merge(sites_dedup, on="sequence.id", how="inner")
    print(f"Peptides with annotations: {len(heatmap_peptide)}")

    # Step 2: Merge peptides at same phosphosite by median
    # Group by gene + phospho.position + module (= unique phosphosite within module)
    grouped = heatmap_peptide.groupby(
        ["Gene.short", "phospho.position", "primary_module", "name_type"],
        dropna=False,
    )

    heatmap_data = grouped.agg(
        n_peptides_merged=("sequence.id", "count"),
        log2FC_3=("log2FC_3", "median"),
        log2FC_15=("log2FC_15", "median"),
        log2FC_30=("log2FC_30", "median"),
    ).reset_index()

    n_merged = (heatmap_data["n_peptides_merged"] > 1).sum()
    print(f"After merging: {len(heatmap_data)} unique phosphosites")
    print(f"Sites where multiple peptides were merged: {n_merged}")

    # Create row labels -- always show phosphosite
    def make_label(row):
        pp = row["phospho.position"]
        if pd.notna(pp) and str(pp).strip() != "":
            return f"{row['Gene.short']} ({pp})"
        return row["Gene.short"]

    heatmap_data["row_label"] = heatmap_data.apply(make_label, axis=1)

    # Order by 15 min log2FC (highest to lowest) within each module
    heatmap_data = (
        heatmap_data
        .sort_values(
            ["primary_module", "log2FC_15"],
            ascending=[True, False],
        )
        .reset_index(drop=True)
    )

    print("\nPhosphosites per module in final heatmap data:")
    print(heatmap_data.groupby("primary_module").size().sort_values(ascending=False).to_string())
    print()

    return heatmap_data


# ============================================================================
# Heatmap generation
# ============================================================================

def generate_colorbar(max_abs, output_dir):
    """Generate a standalone colorbar SVG for the heatmap color scale."""
    
    import matplotlib.colors as mcolors
    from matplotlib.cm import ScalarMappable
    
    # Create figure for colorbar only
    fig, ax = plt.subplots(figsize=(0.3, 2))
    
    # Create colormap and normalizer matching the heatmaps
    cmap = plt.cm.RdBu_r
    norm = mcolors.Normalize(vmin=-max_abs, vmax=max_abs)
    sm = ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    
    # Create colorbar
    cbar = fig.colorbar(sm, cax=ax, orientation='vertical')
    cbar.set_ticks([-max_abs, 0, max_abs])
    cbar.set_ticklabels([f'-{max_abs}', '0', f'{max_abs}'])
    cbar.ax.tick_params(labelsize=FONT_SIZE)
    cbar.outline.set_linewidth(0.3)
    
    # Label
    cbar.set_label(r'$\log_2$(FC)', fontsize=FONT_SIZE)
    
    fig.patch.set_facecolor("white")
    plt.tight_layout(pad=0)
    
    # Save
    colorbar_path = Path(output_dir) / "colorbar.svg"
    fig.savefig(colorbar_path, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    
    print(f"Saved colorbar: {colorbar_path}")
    return colorbar_path


def generate_module_heatmap(mod_data, mod_name, max_abs, output_dir):
    """Generate a heatmap for one module. Save PDF, PNG, SVG."""

    n_sites = len(mod_data)
    if n_sites < 2:
        print(f"  Skipping {mod_name} -- too few phosphosites ({n_sites})")
        return None

    # Build the matrix
    mat = mod_data[["log2FC_3", "log2FC_15", "log2FC_30"]].copy()
    mat = mat.fillna(0)
    mat = mat.clip(-max_abs, max_abs)  # cap at symmetric limits
    mat.index = mod_data["row_label"].values
    mat.columns = ["3 min", "15 min", "30 min"]

    # ---- Figure dimensions with fixed cell size ----
    n_cols = 3
    heatmap_width_in = n_cols * CELL_SIZE_INCHES
    heatmap_height_in = n_sites * CELL_SIZE_INCHES

    # Fixed margins (in inches) for labels
    left_margin = 0    # space for rotated title on left
    right_margin = 0    # space for gene labels on right
    top_margin = 0     # minimal top margin
    bottom_margin = 0  # minimal bottom margin

    # Total figure size
    page_width = left_margin + heatmap_width_in + right_margin
    page_height = top_margin + heatmap_height_in + bottom_margin

    # ---- Create figure and position axes exactly ----
    fig = plt.figure(figsize=(page_width, page_height))
    
    # Calculate axes position in figure coordinates (0-1)
    ax_left = left_margin / page_width
    ax_bottom = bottom_margin / page_height
    ax_width = heatmap_width_in / page_width
    ax_height = heatmap_height_in / page_height
    
    ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])

    # RdBu reversed (matches R rev(brewer.pal(11, "RdBu")))
    cmap = plt.cm.RdBu_r

    # Line width for borders
    inner_border_width = 0.1
    outer_border_width = inner_border_width

    sns.heatmap(
        mat,
        ax=ax,
        cmap=cmap,
        vmin=-max_abs,
        vmax=max_abs,
        center=0,
        linewidths=inner_border_width,
        linecolor="black",
        cbar=False,  # no colorbar
        square=True,
        xticklabels=False,  # remove x-axis labels (will add in Inkscape)
        yticklabels=False,  # we'll add labels on right side manually
        annot=False,
    )

    # Remove tick marks
    ax.tick_params(axis='both', which='both', length=0)

    # Hide default spines
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Draw a rectangle around the heatmap with correct line width
    from matplotlib.patches import Rectangle
    rect = Rectangle((0, 0), n_cols, n_sites, 
                      fill=False, 
                      edgecolor='black', 
                      linewidth=outer_border_width,
                      clip_on=False)
    ax.add_patch(rect)

    # Style adjustments - no axis labels
    ax.set_ylabel("")
    ax.set_xlabel("")

    # Add gene names on the right side of the heatmap
    for i, label in enumerate(mod_data["row_label"].values):
        ax.text(n_cols + 0.15, i + 0.5, label, 
                ha='left', va='center', fontsize=FONT_SIZE_LABELS)

    # Add title on the left side (rotated vertically) with gene count
    title_text = f"{mod_name} ({n_sites})"
    ax.text(-0.2, n_sites / 2, title_text, 
            ha='right', va='center', fontsize=FONT_SIZE,
            rotation=90, transform=ax.transData)

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    
    # DO NOT use tight_layout - it would rescale the axes

    # ---- Save with fixed figure size (no bbox_inches="tight") ----
    base_filename = f"{slugify(mod_name)}_heatmap"

    fig.savefig(Path(output_dir) / f"{base_filename}.pdf",
                dpi=DPI, facecolor="white")
    fig.savefig(Path(output_dir) / f"{base_filename}.png",
                dpi=DPI, facecolor="white")
    fig.savefig(Path(output_dir) / f"{base_filename}.svg",
                facecolor="white")

    plt.close(fig)

    return base_filename


# ============================================================================
# Export summary
# ============================================================================

def export_summary(heatmap_data, output_dir):
    """Export CSV summary of all phosphosites in the heatmap data."""

    export_data = (
        heatmap_data
        .sort_values(["primary_module", "log2FC_15"], ascending=[True, False])
        [["Gene.short", "phospho.position", "row_label", "primary_module",
          "n_peptides_merged", "log2FC_3", "log2FC_15", "log2FC_30"]]
    )

    csv_path = Path(output_dir) / "figure_4D_phosphosite_summary.csv"
    export_data.to_csv(csv_path, index=False)
    print(f"\nExported phosphosite summary: {csv_path}")

    # Module-level summary
    module_summary = (
        heatmap_data
        .groupby("primary_module")
        .agg(
            n_phosphosites=("row_label", "count"),
            n_genes=("Gene.short", "nunique"),
            mean_log2FC_3=("log2FC_3", "mean"),
            mean_log2FC_15=("log2FC_15", "mean"),
            mean_log2FC_30=("log2FC_30", "mean"),
        )
        .sort_values("n_phosphosites", ascending=False)
        .reset_index()
    )

    print("\nFigure 4D module summary:")
    print(module_summary.to_string(index=False))

    return csv_path


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 60)
    print("Figure 4D -- Phosphoproteomic Module Heatmaps")
    print("=" * 60)
    print()

    # Load data
    assignments, res, mdata = load_data()

    # Identify significant phosphosites
    significant_sites = identify_significant_sites(res)

    # Build gene mapping
    gene_mapping = build_gene_mapping(res)

    # Filter to Figure 4D genes + significant sites
    sites_fig4d = filter_figure_4d_sites(gene_mapping, significant_sites, assignments)

    # Compute log2(ctrl.ratio) matrix
    log2fc_wide = compute_log2fc_matrix(mdata, sites_fig4d)

    # Merge peptides at same phosphosite and create row labels
    heatmap_data = merge_and_label(log2fc_wide, sites_fig4d)

    # Use fixed symmetric color scale
    max_abs = MAX_ABS_CAP
    print(f"Color scale: symmetric +/- {max_abs:.2f}")

    # Generate standalone colorbar
    generate_colorbar(max_abs, dir_out)

    # Generate per-module heatmaps
    modules_list = sorted(heatmap_data["primary_module"].unique())
    print(f"\nGenerating heatmaps for {len(modules_list)} Figure 4D modules\n")

    for mod in modules_list:
        mod_data = heatmap_data[heatmap_data["primary_module"] == mod].copy()
        mod_data = mod_data.reset_index(drop=True)
        n_sites = len(mod_data)
        print(f"Module: {mod} ({n_sites} phosphosites)")

        saved = generate_module_heatmap(
            mod_data=mod_data,
            mod_name=mod,
            max_abs=max_abs,
            output_dir=dir_out,
        )
        if saved is not None:
            print(f"  Saved: {saved} (.pdf/.svg/.png)")
        print()

    # Export summary CSV
    export_summary(heatmap_data, dir_out)

    # Final summary
    print()
    print("=" * 60)
    print("Figure 4D Heatmaps Summary")
    print("=" * 60)
    print(f"Sample type: Raw phospho")
    print(f"Gene selection: Curated Figure 4D genes with primary module assignment")
    print(f"Total phosphosites: {len(heatmap_data)}")
    print(f"Unique genes: {heatmap_data['Gene.short'].nunique()}")
    print(f"Primary modules: {heatmap_data['primary_module'].nunique()}")
    print(f"Output directory: {dir_out}")
    print(f"Formats: PDF, SVG, PNG")


if __name__ == "__main__":
    main()