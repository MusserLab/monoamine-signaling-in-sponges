"""
Scatter plot comparing phosphosite-level log2FC values between
tryptamine (15 min, raw phospho) and NO phosphoproteomics.

Generates module-colored dot scatter plot as SVG:
  Dots colored by module (NO signaling, Actin dynamics, Rac1 GTPase signaling)
  for phosphosites significant in both datasets.

Required files (relative to repo root):
  1. phosphoproteomics/data/NO_phosphosites_mapped.csv
  2. phosphoproteomics/data/limma_results_annotated.tsv
  3. phosphoproteomics/data/modules_long_RAW_PHOSPHO_RZ_merged.tsv
  4. data/spongilla_gene_names_final.tsv
  5. phosphoproteomics/data/gene_list_with_modules_Fig4D_Fig4E.tsv

Usage:
  conda env create -f environment.yml && conda activate monoamine-sponges
  python phosphoproteomics/NO_tryptamine_scatter.py
"""

import math
import re

import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from scipy import stats

# ============================================================================
# Configuration
# ============================================================================

W = 2.671  # inches
H = 2.141  # inches
DPI = 300

# Font — use Arial, fall back to Helvetica or sans-serif
for font_name in ["Arial", "Helvetica"]:
    if any(font_name in f.name for f in fm.fontManager.ttflist):
        plt.rcParams["font.family"] = font_name
        break
else:
    plt.rcParams["font.family"] = "sans-serif"

plt.rcParams.update({
    'font.family': 'Arial',
    "font.size": 6,
    "xtick.labelsize": 5,
    "ytick.labelsize": 5,
    "svg.fonttype": "none",  # keep text as text in SVG (not paths)
    'xtick.major.width': 0.5,
    'ytick.major.width': 0.5
})

# Colors (matches R custom_cols4)
COLORS4 = ["#D81B60", "#1E88E5", "#FFC107", "#004D40"]

# Target modules for module-colored dots
MODULES_C = ["NO signaling", "Actin dynamics", "Rac1 GTPase signaling"]
MODULE_COLORS_C = {
    "NO signaling": COLORS4[1],
    "Actin dynamics": COLORS4[2],
    "Rac1 GTPase signaling": COLORS4[3],
    "Other": "grey",
}

# ============================================================================
# Paths
# ============================================================================

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent

dir_data = _SCRIPT_DIR / "data"
dir_annot = _REPO_ROOT / "data"
dir_no = _SCRIPT_DIR / "data"
dir_out = _SCRIPT_DIR / "outs" / "NO_tryptamine_scatter"

dir_out.mkdir(parents=True, exist_ok=True)


# ============================================================================
# Utility functions
# ============================================================================

def normalize_name(x):
    """Remove unicode apostrophes for join matching."""
    if pd.isna(x):
        return x
    return str(x).replace("\u2019", "").replace("'", "")


def slugify(text):
    """Convert text to filename-safe slug."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


def save_figure(fig, filepath_stem):
    """Save figure as SVG."""
    fig.savefig(f"{filepath_stem}.svg", bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Data loading
# ============================================================================

def load_no_data():
    """Load NO phosphosite data with deduplication."""
    path = dir_no / "NO_phosphosites_mapped.csv"

    df = pd.read_csv(path)
    print(f"NO phosphosites loaded: {len(df)} rows")

    # Create site_id and filter out missing positions
    df["site_id"] = df["Protein.ID"].astype(str) + "_" + df["phospho_position"].astype(str)
    df = df[df["phospho_position"].notna() & (df["phospho_position"].astype(str) != "")]

    # Deduplicate: keep the entry with best hit rank per site_id
    hit_rank_map = {"hit": 1, "candidate": 2}
    df["hit_rank"] = df["hit_annotation"].map(hit_rank_map).fillna(3).astype(int)
    df = df.sort_values("hit_rank").drop_duplicates(subset="site_id", keep="first")

    no = df[["site_id", "Protein.ID", "phospho_position",
             "annotation", "logFC", "fdr.limma", "hit_annotation"]].copy()
    no = no.rename(columns={
        "annotation": "Gene_NO",
        "logFC": "logFC_NO",
        "fdr.limma": "fdr_NO",
        "hit_annotation": "hit_annotation_NO",
    })

    print(f"NO unique sites: {len(no)}")
    print(f"NO hits: {(no['hit_annotation_NO'] == 'hit').sum()}")
    print(f"NO candidates: {(no['hit_annotation_NO'] == 'candidate').sum()}")
    return no


def load_tryptamine_data():
    """Load tryptamine raw phospho limma results for 15 min comparison."""
    path_res = dir_data / "limma_results_annotated.tsv"

    res = pd.read_csv(path_res, sep="\t", low_memory=False)

    # Update gene names from final annotations
    path_annot = dir_annot / "spongilla_gene_names_final.tsv"

    annot = pd.read_csv(path_annot, sep="\t",
                        usecols=["Trinity_geneID", "Zang_et_al_2026", "name_type"])
    annot = annot.rename(columns={
        "Trinity_geneID": "trinity_gene_id",
        "Zang_et_al_2026": "Gene_short_new",
        "name_type": "name_type_new",
    })
    annot = annot.drop_duplicates(subset="trinity_gene_id")

    res = res.merge(annot, on="trinity_gene_id", how="left")
    res["Gene.short"] = res["Gene_short_new"].fillna(res["Gene.short"])
    res["name_type"] = res["name_type_new"].fillna("trinity_only")
    res.drop(columns=["Gene_short_new", "name_type_new"], inplace=True)

    # Filter: raw phospho, 15 min comparison
    mask = ((res["sample"] == "phospho") &
            (res["comparison"] == "Tryptamine_15min - DMSO_15min"))
    tryp = res.loc[mask].copy()

    # Create site_id from Protein.ID + phospho.position
    tryp["site_id"] = (tryp["Protein.ID"].astype(str) + "_" +
                       tryp["phospho.position"].astype(str))

    tryp = tryp[["site_id", "sequence.id", "Gene", "Gene.short",
                 "automated_name", "logFC", "fdr", "hit_annotation"]].copy()
    tryp = tryp.rename(columns={
        "Gene": "Gene_Tryp",
        "logFC": "logFC_Tryp",
        "fdr": "fdr_Tryp",
        "hit_annotation": "hit_annotation_Tryp",
    })

    print(f"\nTryptamine raw phospho 15min: {len(tryp)} rows")
    print(f"Unique sites: {tryp['site_id'].nunique()}")
    return tryp


def load_fig4_genes():
    """Load gene list with module boolean columns for Version C."""
    path = dir_data / "gene_list_with_modules_Fig4D_Fig4E.tsv"
    return pd.read_csv(path, sep="\t")


# ============================================================================
# Data joining & classification
# ============================================================================

def join_datasets(no_data, tryp_data):
    """Inner join NO and tryptamine data on site_id; classify significance."""
    merged = no_data.merge(tryp_data, on="site_id", how="inner")
    print(f"Shared phosphosites (detected in both): {len(merged)}")

    merged["sig_NO"] = merged["hit_annotation_NO"].isin(["hit", "candidate"])
    merged["sig_Tryp"] = merged["hit_annotation_Tryp"].isin(["hit", "candidate"])
    merged["sig_either"] = merged["sig_NO"] | merged["sig_Tryp"]
    merged["sig_both"] = merged["sig_NO"] & merged["sig_Tryp"]

    conditions = [
        merged["sig_both"],
        merged["sig_NO"] & ~merged["sig_Tryp"],
        ~merged["sig_NO"] & merged["sig_Tryp"],
    ]
    choices = ["Significant in both", "NO only", "Tryptamine only"]
    merged["sig_category"] = np.select(conditions, choices, default="Neither")

    print("\nSignificance breakdown:")
    print(merged["sig_category"].value_counts().to_string())
    return merged


def compute_stats(plot_data):
    """Compute Pearson correlation and linear regression."""
    x = plot_data["logFC_Tryp"].values
    y = plot_data["logFC_NO"].values

    r, p = stats.pearsonr(x, y)
    # Format with superscript for p-value using mathtext
    if p < 2.2e-16:
        cor_label = r"r = {:.3f}".format(r) + "\n" + r"p < 2.2 $\times$ 10$^{-16}$"
    elif p < 0.001:
        # Extract exponent for superscript formatting
        exp = int(np.floor(np.log10(p)))
        cor_label = f"r = {r:.3f}, p < 10$^{{{exp}}}$"
    else:
        cor_label = f"r = {r:.3f}, p = {p:.3f}"
    print(f"Correlation: r = {r:.3f}, p = {p:.2e}")

    # Axis limits: y symmetric, x starts at -2
    axis_max_raw = max(abs(x).max(), abs(y).max())
    axis_max = math.ceil(axis_max_raw * 1.1)
    print(f"Symmetric axis limit: +/- {axis_max}")

    # Linear regression line from -2 to axis_max
    slope, intercept = np.polyfit(x, y, 1)
    lm_x = np.array([-2, axis_max])
    lm_y = slope * lm_x + intercept

    return cor_label, axis_max, lm_x, lm_y


# ============================================================================
# Shared scatter setup
# ============================================================================

def _setup_scatter_axes(ax, axis_max, cor_label, lm_x, lm_y,
                        title, subtitle):
    """Apply shared formatting: reference lines, regression, labels, axes."""
    # Reference lines
    ax.axhline(0, linestyle="--", color="grey", linewidth=0.3, zorder=0)
    ax.axvline(0, linestyle="--", color="grey", linewidth=0.3, zorder=0)

    # Linear regression line (linewidth 0.5 pt)
    ax.plot(lm_x, lm_y, color="#4D4D4D", linewidth=0.5, zorder=1)

    # Correlation annotation (upper-left) with proper formatting
    ax.text(-1.9, axis_max * 0.95, cor_label,
            ha="left", va="top", fontsize=6, color="black")

    # Axis settings — x starts at -2, y symmetric
    ax.set_xlim(-2, axis_max)
    ax.set_ylim(-axis_max, axis_max)

    # Axis labels with subscript for log2
    ax.set_xlabel(r"Tryptamine 15 min $\log_2$(FC)", fontsize=6)
    ax.set_ylabel(r"NOC-12 3 min $\log_2$(FC)", fontsize=6)

    # No title (will be added in Inkscape if needed)
    # ax.set_title(title, fontsize=7, ha="center")

    # Enclose plotting area with box (all 4 spines visible, linewidth 0.5 pt)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(0.5)
        spine.set_color("black")

    ax.grid(False)
    ax.set_facecolor("white")


# ============================================================================
# Module-colored dot scatter
# ============================================================================

def _get_module_assignments_c(fig4_genes, filter_col=None):
    """Assign priority-based module from fig4 gene list, optionally filtered."""
    gdf = fig4_genes.copy()
    if filter_col is not None:
        gdf = gdf[gdf[filter_col].isin([True, "TRUE", "True", "T", 1, "1"])]

    # Pivot module boolean columns to long format
    rows = []
    for mod in MODULES_C:
        if mod not in gdf.columns:
            continue
        sub = gdf[gdf[mod].isin([True, "TRUE", "True", "T", 1, "1"])][["automated_name"]].copy()
        sub["module_c"] = mod
        sub["priority"] = MODULES_C.index(mod)
        rows.append(sub)

    if not rows:
        return pd.DataFrame(columns=["automated_name", "module_c", "automated_name_norm"])

    long = pd.concat(rows, ignore_index=True)
    # Keep highest priority module per gene
    long = long.sort_values("priority").drop_duplicates(subset="automated_name", keep="first")
    long["automated_name_norm"] = long["automated_name"].apply(normalize_name)
    return long[["automated_name", "module_c", "automated_name_norm"]]


def plot_version_c(plot_data, axis_max, cor_label, lm_x, lm_y,
                   mod_assign, subtitle_suffix):
    """Scatter with module-colored dots for sig-in-both genes."""
    fig, ax = plt.subplots(figsize=(W, H))
    fig.patch.set_facecolor("white")

    pd_copy = plot_data.copy()
    pd_copy["automated_name_norm"] = pd_copy["automated_name"].apply(normalize_name)
    pd_copy = pd_copy.merge(
        mod_assign[["automated_name_norm", "module_c"]],
        on="automated_name_norm", how="left"
    )
    # Only color module genes that are significant in both
    pd_copy["module_c"] = np.where(
        pd_copy["module_c"].notna() & pd_copy["sig_both"],
        pd_copy["module_c"], "Other"
    )

    n_mod = (pd_copy["module_c"] != "Other").sum()

    # Layer 1: Not significant in both (light grey) - increased size
    ns = pd_copy[~pd_copy["sig_both"]]
    ax.scatter(ns["logFC_Tryp"], ns["logFC_NO"],
               c="#CCCCCC", s=1.5, alpha=0.3, linewidths=0,
               zorder=1, rasterized=True)

    # Layer 2: Significant in both but not in a module (dark grey) - increased size
    both_other = pd_copy[pd_copy["sig_both"] & (pd_copy["module_c"] == "Other")]
    ax.scatter(both_other["logFC_Tryp"], both_other["logFC_NO"],
               c="#737373", s=3, alpha=0.5, linewidths=0,
               zorder=2, rasterized=True)

    # Layer 3: Module-colored (significant in both) - increased size
    from matplotlib.lines import Line2D
    legend_handles = []
    for mod_name in MODULES_C:
        mod_pts = pd_copy[pd_copy["module_c"] == mod_name]
        if len(mod_pts) == 0:
            continue
        color = MODULE_COLORS_C[mod_name]
        ax.scatter(mod_pts["logFC_Tryp"], mod_pts["logFC_NO"],
                   c=color, s=6, alpha=0.9, linewidths=0,
                   zorder=3)
        legend_handles.append(
            Line2D([0], [0], marker="o", color="w", linestyle='',
                   markerfacecolor=color, markersize=3, label=mod_name)
        )

    title = "Tryptamine (15 min) vs NO phosphoproteomics"
    subtitle = (f"n = {len(pd_copy)} phosphosites; "
                f"{n_mod} colored ({subtitle_suffix})")
    _setup_scatter_axes(ax, axis_max, cor_label, lm_x, lm_y, title, subtitle)

    # Compact legend in 4th quadrant (lower right, where data is sparse)
    if legend_handles:
        leg = ax.legend(handles=legend_handles, loc="lower right",
                  bbox_to_anchor=(0.98, 0.02),
                  frameon=True, framealpha=1.0, edgecolor='none',
                  fontsize=6, title="Module", title_fontproperties={'weight': 'bold', 'size': 5},
                  handletextpad=0.1, labelspacing=0.2, borderpad=0.3)
        leg.get_frame().set_facecolor('white')

    plt.tight_layout()
    return fig


# ============================================================================
# Summary statistics
# ============================================================================

def print_summary(plot_data):
    """Print summary statistics matching the R version."""
    print("\n=== NO vs Tryptamine 15min Scatter Summary ===\n")
    print(f"Total shared phosphosites: {len(plot_data)}")
    print(f"Significant in either: {plot_data['sig_either'].sum()}")
    print(f"Significant in both: {plot_data['sig_both'].sum()}")

    shared_sig = plot_data[plot_data["sig_both"]].copy()
    if len(shared_sig) > 0:
        print("\nDirection concordance (shared significant sites):")
        both_up = ((shared_sig["logFC_NO"] > 0) & (shared_sig["logFC_Tryp"] > 0)).sum()
        both_down = ((shared_sig["logFC_NO"] < 0) & (shared_sig["logFC_Tryp"] < 0)).sum()
        opposite = len(shared_sig) - both_up - both_down
        print(f"  Both up:   {both_up}")
        print(f"  Both down: {both_down}")
        print(f"  Opposite:  {opposite}")
        pct_same = (both_up + both_down) / len(shared_sig) * 100
        print(f"\nSame direction: {pct_same:.1f}%")

    x = plot_data["logFC_Tryp"].values
    y = plot_data["logFC_NO"].values
    r, p = stats.pearsonr(x, y)
    print(f"\nPearson correlation: r = {r:.3f}, p = {p:.3g}")
    print(f"\nOutput directory: {dir_out}")


# ============================================================================
# Main
# ============================================================================

def main():
    print("Loading data...")

    # 1. Load datasets
    no_data = load_no_data()
    tryp_data = load_tryptamine_data()
    fig4_genes = load_fig4_genes()

    # 2. Join and classify
    plot_data = join_datasets(no_data, tryp_data)

    # 3. Compute statistics
    cor_label, axis_max, lm_x, lm_y = compute_stats(plot_data)

    # ------- Module-colored dots scatter -------
    print("\n=== Module-colored dots ===")
    ma = _get_module_assignments_c(fig4_genes, filter_col=None)
    if len(ma) > 0:
        print(ma["module_c"].value_counts().to_string())

    fig = plot_version_c(plot_data, axis_max, cor_label, lm_x, lm_y,
                         ma, "all module genes")
    save_figure(fig, str(dir_out / "NO_vs_tryptamine_15min_scatter_module_dots"))
    print("Saved (SVG)")

    # ------- Summary -------
    print_summary(plot_data)

    print("\nDone!")
    print(f"Output: {dir_out}")


if __name__ == "__main__":
    main()