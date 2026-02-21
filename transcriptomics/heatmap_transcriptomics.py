"""
Per-module heatmaps for tryptamine transcriptomics.

Generates two versions of per-module heatmaps as SVG:
  1. *_all — heatmap with all genes in the module
  2. *_symbols_only — heatmap restricted to genes with symbol names

Uses DESeq2 Log2FoldChange values directly from the master differential
expression results (not recomputed from raw counts).

Required files (relative to repo root):
  1. transcriptomics/data/Master_DE_Results_Interaction_V2.tsv
  2. transcriptomics/data/transcriptomics_module_assignments.csv
  3. data/spongilla_gene_names_final.tsv

Usage:
  conda env create -f environment.yml && conda activate monoamine-sponges
  python transcriptomics/heatmap_transcriptomics.py
"""

import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)

# ============================================================================
# Configuration
# ============================================================================

MAX_ABS_CAP = 5.0         # cap symmetric color scale at +/- 5
CELL_SIZE_INCHES = 0.112  # each square is 0.095 inch x 0.095 inch
FONT_SIZE = 7             # title font size
FONT_SIZE_LABELS = 6      # gene name labels font size

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': FONT_SIZE,
    'axes.labelsize': FONT_SIZE,
    'xtick.labelsize': FONT_SIZE,
    'ytick.labelsize': FONT_SIZE,
    'svg.fonttype': 'none',  # keep text as text in SVG (not paths)
    'pdf.fonttype': 42,      # TrueType in PDF
})

# Timepoint column order (chronological)
TIMEPOINT_ORDER = ["1h", "4h", "Recovery"]

DPI = 300

# ============================================================================
# Paths  (relative to this script's location in the repo)
# ============================================================================

SCRIPT_DIR = Path(__file__).resolve().parent          # transcriptomics/
REPO_ROOT  = SCRIPT_DIR.parent                        # monoamine-signaling-in-sponges/

dir_data  = SCRIPT_DIR / "data"                       # transcriptomics/data/
dir_annot = REPO_ROOT  / "data"                       # data/  (shared gene names)
dir_out   = SCRIPT_DIR / "outs" / "heatmaps"

dir_out.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Utilities
# ============================================================================

def slugify(text: str) -> str:
    """Convert text to filename-safe slug."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def capitalize_gene_name(name: str) -> str:
    """Capitalize gene name: first letter uppercase, rest lowercase.
    
    E.g., 'GAPDH' -> 'Gapdh', 'ATP1A1' -> 'Atp1a1'
    Preserves any parenthetical suffixes like ' (c12345-g1)'.
    """
    if not name:
        return name
    
    # Check for parenthetical suffix (e.g., " (c12345-g1)")
    match = re.match(r'^(.+?)(\s*\([^)]+\))$', name)
    if match:
        base_name = match.group(1)
        suffix = match.group(2)
        return base_name.capitalize() + suffix
    else:
        return name.capitalize()


def standardize_module(name: str) -> str:
    """Fix typos and inconsistencies in module names."""
    low = name.lower()
    if "vesicular" in low or "veiscular" in low:
        return "Vesicular trafficking"
    if re.search(r"membrane.*ach", low):
        return "Membrane anchoring"
    if "no signal" in low:
        return "NO signaling"
    return name


# ============================================================================
# Data loading
# ============================================================================

def load_master(path_master, path_annot):
    """Load master expression file, join with gene annotations."""
    master = pd.read_csv(path_master, sep="\t", low_memory=False)
    print(f"Master file: {master.shape[0]} rows x {master.shape[1]} cols")
    print(f"Unique genes: {master['GeneID'].nunique()}")
    print(f"Comparisons: {', '.join(master['Comparison_Label'].unique())}")

    # Gene annotations (final names with paralog suffixes)
    annot = pd.read_csv(
        path_annot, sep="\t",
        usecols=["Trinity_geneID", "Zang_et_al_2026", "name_type"],
    ).rename(columns={
        "Trinity_geneID": "trinity_gene_id",
        "Zang_et_al_2026": "Gene_short_new",
    })
    annot = annot.drop_duplicates(subset="trinity_gene_id")

    # Extract Trinity gene ID prefix from GeneID (strip everything after first space)
    # GeneID format: "c100-g1 ..." — Trinity_geneID in gene names file also uses dashes
    master["trinity_gene_id"] = master["GeneID"].str.split(" ").str[0]

    # Join and create Gene.short from final gene names (fall back to Trinity ID)
    master = master.merge(annot, on="trinity_gene_id", how="left")
    master["Gene.short"] = master["Gene_short_new"].fillna(master["trinity_gene_id"])
    master["name_type"] = master["name_type"].fillna("trinity_only")
    master.drop(columns=["Gene_short_new"], inplace=True)

    print(f"Gene name types: {', '.join(sorted(master['name_type'].unique()))}")
    return master




def load_modules(path_modules):
    """Load and standardize module annotations."""
    modules = pd.read_csv(path_modules)
    print(f"Module file: {modules.shape[0]} rows")
    print(f"Genes with module annotation: {modules['Module'].notna().sum()}")

    # Keep only annotated genes and standardize names
    modules = modules[modules["Module"].notna()].copy()
    modules["Module"] = modules["Module"].apply(standardize_module)

    print("\nModule counts after standardization:")
    for mod, cnt in modules["Module"].value_counts().items():
        print(f"  {mod}: {cnt}")

    return modules[["GeneID", "Gene.short", "Module", "Important"]]


# ============================================================================
# Compute log2 fold changes
# ============================================================================

def extract_deseq2_log2fc(master: pd.DataFrame):
    """Extract DESeq2 Log2FoldChange per gene per timepoint.

    Uses the pre-computed Log2FoldChange from DESeq2 output rather than
    recomputing from raw counts, as DESeq2 applies size factor normalization
    and shrinkage for more robust estimates.

    Returns DataFrame with columns:
        GeneID, Gene.short, name_type, log2FC_1h, log2FC_4h, log2FC_Recovery
    """
    # Map comparison labels to timepoint names
    comp_to_tp = {
        "Tryptamine_vs_DMSO_1h": "1h",
        "Tryptamine_vs_DMSO_4h": "4h",
        "Tryptamine_vs_DMSO_Rec": "Recovery",
    }

    # Filter to the 3 direct Tryptamine vs DMSO comparisons
    mask = master["Comparison_Label"].isin(comp_to_tp.keys())
    df = master.loc[mask, ["GeneID", "Gene.short", "name_type",
                           "Comparison_Label", "Log2FoldChange"]].copy()
    df["timepoint"] = df["Comparison_Label"].map(comp_to_tp)
    df["Log2FoldChange"] = pd.to_numeric(df["Log2FoldChange"], errors="coerce")

    print(f"Rows with DESeq2 log2FC: {df['Log2FoldChange'].notna().sum()} / {len(df)}")

    # Pivot timepoints to columns
    log2fc_wide = df.pivot_table(
        index=["GeneID", "Gene.short", "name_type"],
        columns="timepoint",
        values="Log2FoldChange",
    ).reset_index()
    log2fc_wide.columns.name = None

    # Rename columns
    log2fc_wide = log2fc_wide.rename(columns={
        "1h": "log2FC_1h",
        "4h": "log2FC_4h",
        "Recovery": "log2FC_Recovery",
    })

    print(f"Log2FC matrix: {len(log2fc_wide)} genes x 3 timepoints")
    print(f"  Name type distribution: {dict(log2fc_wide['name_type'].value_counts())}")
    return log2fc_wide


# ============================================================================
# Heatmap generation
# ============================================================================

def generate_colorbar(max_abs, output_dir):
    """Generate a standalone colorbar SVG for the heatmap color scale."""
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
    plt.tight_layout(pad=0.1)
    
    # Save
    colorbar_path = output_dir / "colorbar.svg"
    fig.savefig(str(colorbar_path), bbox_inches="tight", facecolor="white")
    plt.close(fig)
    
    print(f"Saved colorbar: {colorbar_path}")
    return colorbar_path


def generate_per_module_heatmap(
    mod_data: pd.DataFrame,
    mod_name: str,
    max_abs: float,
    output_dir: Path,
    suffix: str = "",
):
    """Generate and save a per-module heatmap (PDF, PNG, SVG).

    Args:
        mod_data: DataFrame with Gene.short, log2FC_1h, log2FC_4h, log2FC_Recovery
        mod_name: Module name for title
        max_abs: Symmetric color scale limit
        output_dir: Directory to save to
        suffix: Filename suffix (e.g., "_all", "_symbols_only")

    Returns:
        base_filename or None if too few genes
    """
    from matplotlib.patches import Rectangle

    n_genes = len(mod_data)
    if n_genes < 2:
        return None

    # Build matrix
    fc_cols = ["log2FC_1h", "log2FC_4h", "log2FC_Recovery"]
    mat = mod_data[fc_cols].values.copy()
    mat = np.nan_to_num(mat, nan=0.0)
    mat = np.clip(mat, -max_abs, max_abs)
    
    row_labels = mod_data["Gene.short"].tolist()

    # Capitalize gene names (first letter uppercase, rest lowercase)
    row_labels = [capitalize_gene_name(label) for label in row_labels]

    # Sort rows: Inflammation by 4h, all others by Recovery (highest to lowest)
    if mod_name == "Inflammation":
        sort_col_idx = 1  # log2FC_4h
    else:
        sort_col_idx = 2  # log2FC_Recovery
    order = np.argsort(mat[:, sort_col_idx])[::-1]  # descending
    mat = mat[order]
    row_labels = [row_labels[i] for i in order]

    # Create DataFrame for seaborn
    mat_df = pd.DataFrame(mat, index=row_labels, columns=["1h", "4h", "Recovery"])

    # ---- Figure dimensions with fixed cell size ----
    n_cols = 3
    heatmap_width_in = n_cols * CELL_SIZE_INCHES
    heatmap_height_in = n_genes * CELL_SIZE_INCHES

    # Fixed margins (in inches) for labels
    left_margin = 0
    right_margin = 0
    top_margin = 0
    bottom_margin = 0

    # Total figure size
    page_width = left_margin + heatmap_width_in + right_margin
    page_height = top_margin + heatmap_height_in + bottom_margin

    # ---- Create figure and position axes exactly ----
    fig = plt.figure(figsize=(page_width, page_height))
    
    # Calculate axes position in figure coordinates (0-1)
    ax_left = left_margin / page_width if page_width > 0 else 0
    ax_bottom = bottom_margin / page_height if page_height > 0 else 0
    ax_width = heatmap_width_in / page_width if page_width > 0 else 1
    ax_height = heatmap_height_in / page_height if page_height > 0 else 1
    
    ax = fig.add_axes([ax_left, ax_bottom, ax_width, ax_height])

    # RdBu reversed colormap
    cmap = plt.cm.RdBu_r

    # Line width for borders
    inner_border_width = 0.1
    outer_border_width = inner_border_width

    sns.heatmap(
        mat_df,
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
    rect = Rectangle((0, 0), n_cols, n_genes,
                      fill=False,
                      edgecolor='black',
                      linewidth=outer_border_width,
                      clip_on=False)
    ax.add_patch(rect)

    # Style adjustments - no axis labels
    ax.set_ylabel("")
    ax.set_xlabel("")

    # Add gene names on the right side of the heatmap
    for i, label in enumerate(row_labels):
        ax.text(n_cols + 0.15, i + 0.5, label,
                ha='left', va='center', fontsize=FONT_SIZE_LABELS,
                fontstyle='italic')

    # Add title on the left side (rotated vertically) with gene count
    title_text = f"{mod_name} ({n_genes})"
    ax.text(-0.2, n_genes / 2, title_text,
            ha='right', va='center', fontsize=FONT_SIZE,
            rotation=90, transform=ax.transData)

    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")

    # DO NOT use tight_layout - it would rescale the axes

    # ---- Save with fixed figure size (no bbox_inches="tight") ----
    base_filename = f"{slugify(mod_name)}_heatmap{suffix}"

    fig.savefig(str(output_dir / f"{base_filename}.svg"),
                facecolor="white")

    plt.close(fig)

    return base_filename


# ============================================================================
# Export summaries
# ============================================================================

def export_summary(heatmap_data: pd.DataFrame, output_dir: Path):
    """Export gene-module summary and module-level statistics as CSVs."""
    # Create row labels for export (disambiguate duplicates)
    dup_counts = heatmap_data["Gene.short"].value_counts()
    dups = set(dup_counts[dup_counts > 1].index)

    row_labels = []
    for _, row in heatmap_data.iterrows():
        gs = row["Gene.short"]
        if gs in dups:
            trinity = re.search(r"c\d+-g\d+", row["GeneID"])
            trinity_str = f" ({trinity.group()})" if trinity else ""
            row_labels.append(f"{gs}{trinity_str}")
        else:
            row_labels.append(gs)

    # Capitalize gene names (first letter uppercase, rest lowercase)
    row_labels = [capitalize_gene_name(label) for label in row_labels]

    export = heatmap_data[
        ["GeneID", "Gene.short", "Module", "Important",
         "log2FC_1h", "log2FC_4h", "log2FC_Recovery"]
    ].copy()
    export.insert(2, "row_label", row_labels)
    export = export.sort_values(["Module", "Gene.short"])

    csv_path = output_dir / "transcriptomics_gene_summary.csv"
    export.to_csv(str(csv_path), index=False)
    print(f"\nExported gene-module summary: {csv_path}")

    # Module-level statistics
    fc_cols = ["log2FC_1h", "log2FC_4h", "log2FC_Recovery"]
    mod_stats = (
        heatmap_data
        .groupby("Module")[fc_cols]
        .agg(["count", "mean"])
    )
    # Flatten multi-level columns
    mod_summary = pd.DataFrame({
        "Module": mod_stats.index,
        "n_genes": mod_stats[("log2FC_1h", "count")].values,
        "mean_log2FC_1h": mod_stats[("log2FC_1h", "mean")].values,
        "mean_log2FC_4h": mod_stats[("log2FC_4h", "mean")].values,
        "mean_log2FC_Recovery": mod_stats[("log2FC_Recovery", "mean")].values,
    }).sort_values("n_genes", ascending=False)

    print("\nModule summary statistics:")
    print(mod_summary.to_string(index=False))

    mod_summary.to_csv(str(output_dir / "module_summary_statistics.csv"), index=False)

    return csv_path


# ============================================================================
# Main
# ============================================================================

def main():
    print("=" * 60)
    print("Transcriptomics -- Module Heatmaps")
    print("=" * 60)
    print()

    # --- Load data ---
    print("Loading data...")

    path_master  = dir_data  / "Master_DE_Results_Interaction_V2.tsv"
    path_annot   = dir_annot / "spongilla_gene_names_final.tsv"
    path_modules = dir_data  / "transcriptomics_module_assignments.csv"

    master = load_master(path_master, path_annot)
    modules = load_modules(path_modules)

    # --- Extract DESeq2 log2FC ---
    print("\n" + "=" * 60)
    print("Extracting DESeq2 log2 fold changes")
    print("=" * 60)

    log2fc_wide = extract_deseq2_log2fc(master)

    # --- Join with modules ---
    # Join on GeneID only (modules has old Gene.short values)
    heatmap_data = log2fc_wide.merge(
        modules[["GeneID", "Module", "Important"]],
        on="GeneID",
        how="inner",
    )
    print(f"\nGenes with module annotations and expression data: {len(heatmap_data)}")
    print(f"  Name type distribution: {dict(heatmap_data['name_type'].value_counts())}")
    print(f"\nGenes per module:")
    for mod, cnt in heatmap_data["Module"].value_counts().items():
        print(f"  {mod}: {cnt}")

    # --- Determine global color scale ---
    fc_cols = ["log2FC_1h", "log2FC_4h", "log2FC_Recovery"]
    all_vals = heatmap_data[fc_cols].values.ravel()
    all_vals = all_vals[~np.isnan(all_vals)]
    max_abs = min(np.max(np.abs(all_vals)), MAX_ABS_CAP)
    print(f"\nColor scale: symmetric +/- {max_abs:.2f}")

    # --- Generate standalone colorbar ---
    generate_colorbar(max_abs, dir_out)

    # --- Generate per-module heatmaps (dual versions) ---
    print("\n" + "=" * 60)
    print("Generating per-module heatmaps")
    print("=" * 60)

    modules_list = sorted(heatmap_data["Module"].unique())

    for mod in modules_list:
        mod_data_all = heatmap_data[heatmap_data["Module"] == mod].copy()
        mod_data_symbols = mod_data_all[mod_data_all["name_type"] == "symbol"].copy()

        n_all = len(mod_data_all)
        n_sym = len(mod_data_symbols)
        print(f"\nModule: {mod} ({n_all} total, {n_sym} with symbols)")

        # Version 1: All genes
        if n_all >= 2:
            saved = generate_per_module_heatmap(
                mod_data_all, mod, max_abs, dir_out, suffix="_all",
            )
            if saved:
                print(f"  Saved: {saved} (.svg)")
        else:
            print("  Skipping all version - too few genes")

        # Version 2: Symbols only
        if n_sym >= 2:
            saved = generate_per_module_heatmap(
                mod_data_symbols, mod, max_abs, dir_out, suffix="_symbols_only",
            )
            if saved:
                print(f"  Saved: {saved} (.svg)")
        else:
            print("  Skipping symbols_only version - too few symbol genes")

    # --- Export CSV summaries ---
    print("\n" + "=" * 60)
    print("Exporting summaries")
    print("=" * 60)

    export_summary(heatmap_data, dir_out)

    # --- Final summary ---
    print()
    print("=" * 60)
    print("Transcriptomics Heatmaps Summary")
    print("=" * 60)
    print(f"Total genes: {len(heatmap_data)}")
    print(f"Unique modules: {heatmap_data['Module'].nunique()}")
    print(f"Output directory: {dir_out}")
    print(f"Format: SVG")


if __name__ == "__main__":
    main()