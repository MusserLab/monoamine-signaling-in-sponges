"""
Supplemental per-module heatmaps for raw phosphoproteomics.

Replicates the R script phospho_heatmaps_supplemental.qmd in Python.

Generates per-module heatmaps for selected signaling modules, showing
per-timepoint DMSO-normalized log2FC (log2(Tryptamine_Xmin / DMSO_Xmin))
for all hit/candidate phosphosites. Phosphopeptides mapping to the same
phosphosite are merged by median. Rows ordered by log2FC at 15 min
(highest to lowest).

Modules (in priority order for multi-membership resolution):
  Synaptic (merged from Post-synaptic scaffold + Exocytosis),
  MTORC2 signaling, GPCR signaling, Actomyosin contractility,
  Phospholipid signaling, Rho signaling, Regulation of adhesion complexes

Only genes with Plot_Module == TRUE are included.

Required files (relative to repo root):
  1. phosphoproteomics/data/limma_results_annotated.tsv
  2. phosphoproteomics/data/mdata.rds  (read via pyreadr)
  3. phosphoproteomics/data/modules_long_RAW_PHOSPHO_RZ_merged.tsv
  4. data/spongilla_gene_names_final.tsv

Usage:
  conda env create -f environment.yml && conda activate monoamine-sponges
  python phosphoproteomics/heatmap_supplemental.py
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
CELL_SIZE_INCHES = 0.07  # each square is 0.07 inch x 0.07 inch
FONT_SIZE = 7            # title font size
FONT_SIZE_LABELS = 5     # protein name labels font size
MAX_ROWS_PER_COLUMN = 50 # maximum rows per column before splitting into multiple columns

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE,
    'ytick.labelsize': FONT_SIZE,
    'svg.fonttype': 'none',  # keep text as text in SVG (not paths)
})

# Modules to generate supplemental heatmaps for, in priority order (highest first).
# When a gene belongs to multiple modules, it is assigned to the highest-priority one.
SUPPLEMENTAL_MODULES = [
    "Synaptic",
    "MTORC2 signaling",
    "GPCR signaling",
    "Actomyosin contractility",
    "Phospholipid signaling",
    "Rho signaling",
    "Regulation of adhesion complexes",
]

# The 5 key comparisons (using the "comparison" column format)
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
dir_out = _REPO_ROOT / "outs" / "phosphoproteomics" / "supplemental_heatmaps"

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
    """Load limma results, gene annotations, modules, and mdata."""
    print("Loading limma results...")
    path_res = dir_data / "limma_results_annotated.tsv"
    res = pd.read_csv(path_res, sep="\t", low_memory=False)

    # Gene annotations — join on automated_name (= seurat_name) for robust matching
    print("Loading gene annotations...")
    path_annot = dir_annot / "spongilla_gene_names_final.tsv"
    annot = pd.read_csv(
        path_annot,
        sep="\t",
        usecols=["seurat_name", "Zang_et_al_2026", "name_type"],
    ).rename(columns={
        "Zang_et_al_2026": "Gene_short_new",
        "name_type": "name_type_new",
    })
    annot["automated_name_norm"] = annot["seurat_name"].apply(normalize_name)
    annot = annot.drop_duplicates(subset="automated_name_norm")
    annot = annot[["automated_name_norm", "Gene_short_new", "name_type_new"]]

    # Update Gene.short with final gene names (Zang_et_al_2026)
    res["automated_name_norm"] = res["automated_name"].apply(normalize_name)
    res = res.merge(annot, on="automated_name_norm", how="left")
    res["Gene.short"] = res["Gene_short_new"].fillna(res["Gene.short"])
    if "name_type" in res.columns:
        res["name_type"] = res["name_type_new"].fillna(res["name_type"])
    else:
        res["name_type"] = res["name_type_new"].fillna("trinity_only")
    res = res.drop(columns=["Gene_short_new", "name_type_new"])

    print("Loading modules...")
    path_modules = dir_data / "modules_long_RAW_PHOSPHO_RZ_merged.tsv"
    modules_raw = pd.read_csv(path_modules, sep="\t")
    modules = modules_raw[modules_raw["Plot_Module"] == True].copy()
    modules["module"] = modules["module"].replace({
        "Post-synaptic scaffold": "Synaptic",
        "Exocytosis": "Synaptic",
    })
    modules["automated_name_norm"] = modules["Automated.name"].apply(normalize_name)
    modules = modules[["automated_name_norm", "clean_annotation_name_fixed", "module"]].copy()
    modules = modules.rename(columns={"clean_annotation_name_fixed": "gene_name"})
    modules = modules.drop_duplicates()

    print("Loading mdata from RDS (this may take a moment)...")
    path_mdata = dir_data / "mdata.rds"
    mdata = pyreadr.read_r(str(path_mdata))[None]

    print(f"  res: {len(res)} rows")
    print(f"  modules: {len(modules)} rows")
    print(f"  mdata: {len(mdata)} rows")

    return res, modules, mdata


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
    """Generate a heatmap for one module, splitting into multiple columns if needed.
    
    Large modules (>MAX_ROWS_PER_COLUMN) are split into multiple column groups,
    each containing up to MAX_ROWS_PER_COLUMN rows. Each column group has 3 timepoint
    columns (3 min, 15 min, 30 min).
    """
    import math
    from matplotlib.patches import Rectangle

    n_sites = len(mod_data)
    if n_sites < 2:
        print(f"  Skipping {mod_name} -- too few phosphosites ({n_sites})")
        return None

    # Determine number of column groups needed
    n_column_groups = math.ceil(n_sites / MAX_ROWS_PER_COLUMN)
    rows_per_group = math.ceil(n_sites / n_column_groups)  # distribute evenly
    
    # Split data into chunks
    chunks = []
    for i in range(n_column_groups):
        start_idx = i * rows_per_group
        end_idx = min((i + 1) * rows_per_group, n_sites)
        chunk = mod_data.iloc[start_idx:end_idx].copy().reset_index(drop=True)
        chunks.append(chunk)
    
    # Calculate the max rows in any chunk (for consistent height)
    max_rows_in_chunk = max(len(chunk) for chunk in chunks)
    
    # ---- Figure dimensions ----
    n_timepoint_cols = 3  # 3 min, 15 min, 30 min
    heatmap_group_width = n_timepoint_cols * CELL_SIZE_INCHES
    heatmap_height = max_rows_in_chunk * CELL_SIZE_INCHES
    
    # Spacing between column groups (in inches)
    column_gap = 0.9  # gap between heatmap groups for labels
    
    # Total width: all heatmap groups + gaps between them
    total_heatmap_width = (n_column_groups * heatmap_group_width + 
                           (n_column_groups - 1) * column_gap)
    
    # Margins
    left_margin = 0.02
    right_margin = 0.02
    top_margin = 0.02
    bottom_margin = 0.02
    
    page_width = left_margin + total_heatmap_width + right_margin
    page_height = top_margin + heatmap_height + bottom_margin

    # ---- Create figure ----
    fig = plt.figure(figsize=(page_width, page_height))
    
    # RdBu reversed colormap
    cmap = plt.cm.RdBu_r
    inner_border_width = 0.1
    outer_border_width = inner_border_width

    # ---- Draw each column group ----
    for group_idx, chunk in enumerate(chunks):
        chunk_n_rows = len(chunk)
        
        # Build matrix for this chunk
        mat = chunk[["log2FC_3", "log2FC_15", "log2FC_30"]].copy()
        mat = mat.fillna(0)
        mat = mat.clip(-max_abs, max_abs)
        mat.index = chunk["row_label"].values
        mat.columns = ["3 min", "15 min", "30 min"]
        
        # Calculate axes position for this column group
        # X position: left_margin + (group_idx * (heatmap_width + gap))
        ax_left_inches = left_margin + group_idx * (heatmap_group_width + column_gap)
        ax_left = ax_left_inches / page_width
        
        # Y position: align to bottom (all groups same height visually, but may have fewer rows)
        # We want the heatmap to be top-aligned, so calculate from top
        chunk_height_inches = chunk_n_rows * CELL_SIZE_INCHES
        ax_bottom_inches = bottom_margin + (heatmap_height - chunk_height_inches)
        ax_bottom = ax_bottom_inches / page_height
        
        ax_width = heatmap_group_width / page_width
        ax_height = chunk_height_inches / page_height
        
        ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])
        
        # Draw heatmap
        sns.heatmap(
            mat,
            ax=ax,
            cmap=cmap,
            vmin=-max_abs,
            vmax=max_abs,
            center=0,
            linewidths=inner_border_width,
            linecolor="black",
            cbar=False,
            square=True,
            xticklabels=False,
            yticklabels=False,
            annot=False,
        )
        
        # Remove tick marks
        ax.tick_params(axis='both', which='both', length=0)
        
        # Hide default spines
        for spine in ax.spines.values():
            spine.set_visible(False)
        
        # Draw rectangle border around this heatmap chunk
        rect = Rectangle((0, 0), n_timepoint_cols, chunk_n_rows,
                          fill=False,
                          edgecolor='black',
                          linewidth=outer_border_width,
                          clip_on=False)
        ax.add_patch(rect)
        
        # No axis labels
        ax.set_ylabel("")
        ax.set_xlabel("")
        
        # Add gene names on the right side of this column group
        for i, label in enumerate(chunk["row_label"].values):
            ax.text(n_timepoint_cols + 0.15, i + 0.5, label,
                    ha='left', va='center', fontsize=FONT_SIZE_LABELS)

    # Add module title on the left side of the first column group (rotated vertically)
    if mod_name == "Synaptic":
        display_name = "Synaptic signaling"
    else:
        display_name = mod_name
    
    title_text = f"{display_name} ({n_sites})"
    
    # Position title to the left of the first heatmap
    # Use figure coordinates for the title
    first_ax = fig.axes[0]
    # Get the left edge of the first axes in figure coordinates
    ax_bbox = first_ax.get_position()
    title_x = ax_bbox.x0 - 0.01  # slightly to the left of axes
    title_y = ax_bbox.y0 + ax_bbox.height / 2  # vertically centered
    
    fig.text(title_x, title_y, title_text,
             ha='right', va='center', fontsize=FONT_SIZE,
             rotation=90, transform=fig.transFigure)

    fig.patch.set_facecolor("white")

    # ---- Save ----
    base_filename = f"{slugify(display_name)}_heatmap"

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
        .sort_values(["module", "log2FC_15"], ascending=[True, False])
        [["Gene.short", "phospho.position", "gene_name", "row_label", "module",
          "n_peptides_merged", "log2FC_3", "log2FC_15", "log2FC_30"]]
    )

    csv_path = Path(output_dir) / "supplemental_phosphosite_summary.csv"
    export_data.to_csv(csv_path, index=False)
    print(f"\nExported phosphosite summary: {csv_path}")

    # Module-level summary
    module_summary = (
        heatmap_data
        .groupby("module")
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

    print("\nSupplemental module summary:")
    print(module_summary.to_string(index=False))

    return csv_path


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 60)
    print("Supplemental -- Phosphoproteomic Module Heatmaps")
    print("=" * 60)
    print()

    res, modules, mdata = load_data()

    # ---- Filter modules to supplemental set ----
    modules_supp = modules[modules["module"].isin(SUPPLEMENTAL_MODULES)].copy()
    found = sorted(set(SUPPLEMENTAL_MODULES) & set(modules_supp["module"].unique()))
    missing = sorted(set(SUPPLEMENTAL_MODULES) - set(modules_supp["module"].unique()))

    print(f"\nAll modules: {modules['module'].nunique()}")
    print(f"Supplemental modules selected: {modules_supp['module'].nunique()}")
    print(f"  Found: {', '.join(found)}")
    if missing:
        print(f"  MISSING: {', '.join(missing)}")

    # Deduplicate: assign each gene to its highest-priority module
    # Priority is defined by position in SUPPLEMENTAL_MODULES (index 0 = highest)
    priority_map = {mod: i for i, mod in enumerate(SUPPLEMENTAL_MODULES)}
    modules_supp["priority"] = modules_supp["module"].map(priority_map)
    modules_supp = (
        modules_supp
        .sort_values("priority")
        .drop_duplicates(subset="automated_name_norm", keep="first")
        .drop(columns="priority")
    )
    print(f"Genes after priority dedup: {len(modules_supp)} (each gene in exactly one module)")

    print("\nGenes per supplemental module:")
    for mod, count in modules_supp.groupby("module").size().sort_values(ascending=False).items():
        print(f"  {mod}: {count}")

    # ---- Identify hits/candidates across 5 key comparisons ----
    sig_mask = (
        (res["sample"] == SAMPLE_TYPE) &
        (res["comparison"].isin(RELEVANT_COMPARISONS)) &
        (res["hit_annotation"].isin(["hit", "candidate"]))
    )
    significant_sites = res.loc[sig_mask, "sequence.id"].unique()
    print(f"\nSignificant phosphosites (raw phospho): {len(significant_sites)}")

    # ---- Map phosphosites to modules ----
    gene_mapping = (
        res[res["sample"] == SAMPLE_TYPE]
        [["sequence.id", "Gene", "Gene.short", "automated_name_norm",
          "name_type", "phospho.position"]]
        .drop_duplicates()
    )

    sites_sig = gene_mapping[gene_mapping["sequence.id"].isin(significant_sites)]
    sites_with_modules = sites_sig.merge(modules_supp, on="automated_name_norm", how="inner")

    print(f"Phosphosites in supplemental modules: {sites_with_modules['sequence.id'].nunique()}")
    print(f"Unique genes: {sites_with_modules['Gene.short'].nunique()}")

    print("\nPhosphosites per supplemental module:")
    for mod, count in sites_with_modules.groupby("module")["sequence.id"].nunique().sort_values(ascending=False).items():
        print(f"  {mod}: {count}")

    # ---- Compute per-timepoint DMSO-normalized log2FC ----
    print("\nComputing per-timepoint DMSO-normalized log2FC...")

    target_ids = set(sites_with_modules["sequence.id"].unique())
    # Pull ctrl.ratio for BOTH Tryptamine AND DMSO conditions
    tc_data = mdata[
        (mdata["sample"] == SAMPLE_TYPE) &
        (mdata["measurement"] == "ctrl.ratio") &
        (mdata["sequence.id"].isin(target_ids))
    ].copy()

    # Extract timepoint and treatment
    tc_data["time_min"] = tc_data["condition"].str.extract(r"(\d+)").astype(int)
    tc_data["treatment"] = np.where(
        tc_data["condition"].str.startswith("Tryptamine"), "Tryptamine", "DMSO"
    )
    tc_data["log2_ctrl_ratio"] = np.log2(tc_data["value"])

    # Remove non-finite values
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
        index="sequence.id", columns="time_min", values="log2FC"
    ).reset_index()
    log2fc_wide.columns.name = None
    log2fc_wide = log2fc_wide.rename(columns={3: "log2FC_3", 15: "log2FC_15", 30: "log2FC_30"})
    print(f"Before merging peptides: {len(log2fc_wide)} peptides x 3 timepoints")

    # ---- Join with gene/module annotations (still at peptide level) ----
    heatmap_peptide = log2fc_wide.merge(
        sites_with_modules[["sequence.id", "Gene.short", "gene_name", "module",
                            "name_type", "phospho.position"]],
        on="sequence.id",
        how="inner",
    )
    print(f"Peptides with annotations: {len(heatmap_peptide)}")

    # ---- Step 2: Merge peptides at same phosphosite by median ----
    group_cols = ["Gene.short", "phospho.position", "module", "gene_name", "name_type"]
    agg_dict = {
        "sequence.id": "count",    # n_peptides_merged
        "log2FC_3": "median",
        "log2FC_15": "median",
        "log2FC_30": "median",
    }
    heatmap_data = (
        heatmap_peptide
        .groupby(group_cols, dropna=False)
        .agg(agg_dict)
        .reset_index()
        .rename(columns={"sequence.id": "n_peptides_merged"})
    )

    n_merged = (heatmap_data["n_peptides_merged"] > 1).sum()
    print(f"After merging: {len(heatmap_data)} unique phosphosites")
    print(f"Sites where multiple peptides were merged: {n_merged}")

    # ---- Create row labels ----
    def make_row_label(row):
        pp = row["phospho.position"]
        gs = row["Gene.short"]
        if pd.notna(pp) and str(pp).strip() != "":
            return f"{gs} ({pp})"
        return str(gs)

    heatmap_data["row_label"] = heatmap_data.apply(make_row_label, axis=1)

    # ---- Order by log2FC_15 descending within each module ----
    heatmap_data = (
        heatmap_data
        .sort_values(["module", "log2FC_15"], ascending=[True, False])
        .reset_index(drop=True)
    )

    print("\nPhosphosites per module in final heatmap data:")
    for mod, count in heatmap_data.groupby("module").size().sort_values(ascending=False).items():
        print(f"  {mod}: {count}")

    # ---- Use fixed symmetric color scale ----
    max_abs = MAX_ABS_CAP
    print(f"\nColor scale: symmetric +/- {max_abs:.2f}")

    # ---- Generate standalone colorbar ----
    generate_colorbar(max_abs, dir_out)

    # ---- Generate per-module heatmaps ----
    modules_list = sorted(heatmap_data["module"].unique())
    print(f"\nGenerating heatmaps for {len(modules_list)} supplemental modules\n")

    for mod in modules_list:
        mod_data = heatmap_data[heatmap_data["module"] == mod].copy()
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

    # ---- Export summary CSV ----
    export_summary(heatmap_data, dir_out)

    # ---- Final summary ----
    print()
    print("=" * 60)
    print("Supplemental Heatmaps Summary")
    print("=" * 60)
    print(f"Sample type: Raw phospho")
    print(f"Gene selection: Plot_Module genes in supplemental modules")
    print(f"Total phosphosites: {len(heatmap_data)}")
    print(f"Unique genes: {heatmap_data['Gene.short'].nunique()}")
    print(f"Modules: {heatmap_data['module'].nunique()}")
    print(f"Output directory: {dir_out}")
    print(f"Formats: PDF, SVG, PNG")


if __name__ == "__main__":
    main()