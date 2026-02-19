"""
Volcano plots for tryptamine phosphoproteomics (raw phospho).

Generates three types of volcano plots as PDF, PNG, and SVG:
  1. volcano_hits_raw_phospho/ — hits and candidates highlighted (single color)
  2. volcano_module_overlay_raw_phospho/ — multiple modules colored on same plot
  3. volcano_per_module_raw_phospho/ — individual module plots with gene labels

Required files (relative to repo root):
  1. phosphoproteomics/data/limma_results_annotated.tsv
  2. phosphoproteomics/data/modules_long_RAW_PHOSPHO_RZ_merged.tsv
  3. data/spongilla_gene_names_final.tsv

Usage:
  pip install pandas matplotlib adjustText
  python phosphoproteomics/volcano_plot.py
"""

import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from pathlib import Path

# ============================================================================
# Configuration
# ============================================================================

W = 2.224   # inches
H = 1.616   # inches
DPI = 1500   # Higher DPI for sharper rasterized elements

plt.rcParams.update({
    # ---------- Fonts ----------
    'font.family': 'Arial',
    'font.size': 5,
    'axes.titlesize': 6,
    'axes.labelsize': 5,
    'xtick.labelsize': 5,
    'ytick.labelsize': 5,
    'legend.fontsize': 5,
    # ---------- Figure ----------
    'figure.figsize': (W, H),
    # ---------- Padding ----------
    'axes.titlepad': 0.25,     # title → axes padding (points)
    'axes.labelpad': 0.25,     # x/y label → axes padding (points)
    'xtick.major.pad': 0.25,   # x tick label padding (points)
    'ytick.major.pad': 0.25,   # y tick label padding (points)
    # ---------- Tick marks ----------
    'xtick.major.size': 2,     # tick length (points)
    'ytick.major.size': 2,
    'xtick.major.width': 0.3,  # tick thickness
    'ytick.major.width': 0.3,
    # ---------- Axes border (spines) ----------
    'axes.linewidth': 0.3,     # border thickness
    # ---------- SVG ----------
    'svg.fonttype': 'none',    # keep text as text in SVG (not paths)
})

# Colors (matches R custom_cols4)
COLORS4 = ["#D81B60", "#1E88E5", "#FFC107", "#004D40"]
COLOR_BG = "#B3B3B3"  # grey70

# Comparisons
COMPARISONS = [
    "Tryptamine_3min vs DMSO_3min",
    "Tryptamine_15min vs DMSO_15min",
    "Tryptamine_30min vs DMSO_30min",
    "(Tryptamine_15min vs DMSO_15min) against (Tryptamine_3min vs DMSO_3min)",
    "(Tryptamine_30min vs DMSO_30min) against (Tryptamine_15min vs DMSO_15min)",
]

SAMPLE = "phospho"

# Module overlay sets (matches R cytoskeleton_adhesion_modules)
OVERLAY_MODULES = [
    "Actomyosin contractility",
    "GPCR signaling",
    "Rho signaling",
    "Regulation of adhesion complexes",
]

# ============================================================================
# Paths
# ============================================================================


_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _SCRIPT_DIR.parent

dir_data = _SCRIPT_DIR / "data"
dir_annot = _REPO_ROOT / "data"

dir_hits = _REPO_ROOT / "outs" / "phosphoproteomics" / "volcano_hits"
dir_overlay = _REPO_ROOT / "outs" / "phosphoproteomics" / "volcano_overlay"
dir_per_module = _REPO_ROOT / "outs" / "phosphoproteomics" / "volcano_per_module"

# ============================================================================
# Data loading
# ============================================================================

def normalize_name(x):
    """Remove apostrophes for join matching."""
    if pd.isna(x):
        return x
    return str(x).replace("\u2019", "").replace("'", "")


def load_data():
    """Load limma results, gene annotations, and modules."""
    # Limma results
    res = pd.read_csv(dir_data / "limma_results_annotated.tsv", sep="\t",
                      low_memory=False)

    # Gene annotations (final gene names with paralog suffixes)
    annot = pd.read_csv(
        dir_annot / "spongilla_gene_names_final.tsv",
        sep="\t",
        usecols=["Trinity_geneID", "Zang_et_al_2026", "name_type"],
    ).rename(columns={"Trinity_geneID": "trinity_gene_id",
                       "Zang_et_al_2026": "Gene_short_new",
                       "name_type": "name_type_new"})
    annot = annot.drop_duplicates(subset="trinity_gene_id")

    res = res.merge(annot, on="trinity_gene_id", how="left")
    res["Gene.short"] = res["Gene_short_new"].fillna(res["Gene.short"])
    res["name_type"] = res["name_type_new"].fillna("trinity_only")
    res = res.drop(columns=["Gene_short_new", "name_type_new"])

    # Normalized name for module joining
    res["Gene_norm"] = res["automated_name"].apply(normalize_name)

    # Modules
    modules = pd.read_csv(dir_data / "modules_long_RAW_PHOSPHO_RZ_merged.tsv", sep="\t")
    modules["Gene_norm"] = modules["Automated.name"].apply(normalize_name)
    modules["Plot_Module"] = modules["Plot_Module"].isin([True, "TRUE", "True", "T", 1, "1"])
    modules["Label_Name"] = modules["Label_Name"].isin([True, "TRUE", "True", "T", 1, "1"])

    return res, modules


def slugify(text):
    """Convert text to filename-safe slug."""
    return re.sub(r"[^a-z0-9]+", "_", text.lower()).strip("_")


# ============================================================================
# Shared plotting helpers
# ============================================================================

def _setup_axes(ax):
    """Configure axes consistently across all plot types."""
    ax.set_ylim(bottom=0)
    ax.set_xlabel(r"$\log_2$(FC)")
    ax.set_ylabel(r"$-\log_{10}$(p-value)")
    # Fix y-axis label position to be consistent regardless of tick label width
    ax.yaxis.set_label_coords(-0.05, 0.5)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    ax.grid(False)


def _plot_background(ax, bg):
    """Plot grey background points (rasterized for performance)."""
    ax.scatter(bg["logFC"], -np.log10(bg["pvalue"]),
               c=COLOR_BG, s=2, alpha=0.35, marker="o",
               linewidths=0, zorder=1, rasterized=True)


def _setup_figure_dpi(fig):
    """Set rasterization DPI for the figure."""
    fig.set_dpi(DPI)


def save_figure(fig, filepath_stem):
    """Save figure as PDF, PNG, and SVG."""
    fig.savefig(f"{filepath_stem}.pdf", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{filepath_stem}.png", dpi=DPI, bbox_inches="tight")
    fig.savefig(f"{filepath_stem}.svg", bbox_inches="tight")
    plt.close(fig)


# ============================================================================
# Plot type 1: Hits/Candidates
# ============================================================================

def plot_volcano_hits(df, comparison):
    """Volcano with hits (one color) and candidates highlighted."""
    bg = df[df["hit_annotation"] == "no hit"]
    fg = df[df["hit_annotation"].isin(["candidate", "hit"])]

    fig, ax = plt.subplots(figsize=(W, H))
    _setup_figure_dpi(fig)
    _plot_background(ax, bg)

    # All foreground in single color (matches R raw phospho script)
    if len(fg) > 0:
        ax.scatter(fg["logFC"], -np.log10(fg["pvalue"]),
                   c=COLORS4[0], s=4, alpha=0.75, marker="o",
                   linewidths=0, zorder=2)

    _setup_axes(ax)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    plt.tight_layout(pad=0.3)
    return fig


# ============================================================================
# Plot type 2: Multi-module overlay
# ============================================================================

def plot_volcano_overlay(df, modules, modules_keep, comparison):
    """Volcano with multiple modules colored on the same plot."""
    # Assign each gene to its first matching module
    mod_sub = modules[modules["module"].isin(modules_keep)].copy()
    mod_sub["rank"] = mod_sub["module"].apply(lambda m: modules_keep.index(m))
    mod_sub = (mod_sub
               .sort_values("rank")
               .drop_duplicates(subset="Gene_norm", keep="first"))
    mod_gene_to_module = dict(zip(mod_sub["Gene_norm"], mod_sub["module"]))

    df = df.copy()
    df["module"] = df["Gene_norm"].map(mod_gene_to_module).fillna("Other")
    bg = df.copy()
    fg = df[(df["module"] != "Other") & df["hit_annotation"].isin(["hit", "candidate"])]

    fig, ax = plt.subplots(figsize=(W, H))
    _setup_figure_dpi(fig)
    _plot_background(ax, bg)

    # Plot each module in its color
    for i, mod_name in enumerate(modules_keep):
        mod_fg = fg[fg["module"] == mod_name]
        if len(mod_fg) == 0:
            continue
        color = COLORS4[i % len(COLORS4)]
        ax.scatter(mod_fg["logFC"], -np.log10(mod_fg["pvalue"]),
                   c=color, s=4, alpha=0.8, marker="o",
                   linewidths=0, zorder=2)

    _setup_axes(ax)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    plt.tight_layout(pad=0.3)
    return fig


# ============================================================================
# Plot type 3: Per-module with gene labels
# ============================================================================

def plot_volcano_per_module(df, modules, module_name, comparison,
                            mod_color, label_genes=True):
    """Volcano for a single module with optional gene labels."""
    df = df.copy()
    mod_genes = set(modules[modules["module"] == module_name]["Gene_norm"].unique())
    df["in_module"] = df["Gene_norm"].isin(mod_genes)
    bg = df.copy()
    fg = df[df["in_module"] & df["hit_annotation"].isin(["hit", "candidate"])]

    fig, ax = plt.subplots(figsize=(W, H))
    _setup_figure_dpi(fig)
    _plot_background(ax, bg)

    # Module points
    if len(fg) > 0:
        ax.scatter(fg["logFC"], -np.log10(fg["pvalue"]),
                   c=mod_color, s=4, alpha=0.8, marker="o",
                   linewidths=0, zorder=2)

    # Gene labels
    if label_genes and len(fg) > 0:
        # One label per gene (lowest p-value)
        lab = (fg
               .sort_values("pvalue")
               .drop_duplicates(subset="Gene_norm", keep="first"))
        texts = []
        for _, row in lab.iterrows():
            t = ax.text(row["logFC"], -np.log10(row["pvalue"]),
                        row["Gene.short"],
                        fontsize=3, ha="center", va="center",
                        zorder=3)
            texts.append(t)
        if texts:
            adjust_text(texts, ax=ax,
                        arrowprops=dict(arrowstyle="-", color="grey", alpha=0.4, lw=0.3),
                        expand=(1.2, 1.4))

    _setup_axes(ax)
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    plt.tight_layout(pad=0.3)
    return fig


# ============================================================================
# Main
# ============================================================================

def main():
    print("Loading data...")
    res, modules = load_data()
    all_modules = sorted(modules["module"].unique())

    # Create output directories
    dir_hits.mkdir(parents=True, exist_ok=True)
    dir_overlay.mkdir(parents=True, exist_ok=True)
    dir_per_module.mkdir(parents=True, exist_ok=True)

    # ------- 1. Hits/Candidates -------
    print("\n=== Hits/Candidates ===")
    for comp in COMPARISONS:
        df = res[(res["comparison.label"] == comp) & (res["sample"] == SAMPLE)]
        if len(df) == 0:
            print(f"  SKIP: {comp}")
            continue
        n_hits = (df["hit_annotation"] == "hit").sum()
        n_cands = (df["hit_annotation"] == "candidate").sum()
        print(f"  PLOT: {comp[:50]}... | n={len(df)} | hits={n_hits} | cands={n_cands}")
        fig = plot_volcano_hits(df, comp)
        save_figure(fig, str(dir_hits / f"volcano_hits_phospho_{slugify(comp)}"))

    # ------- 2. Module overlay -------
    print("\n=== Module overlay ===")
    print(f"  Modules: {', '.join(OVERLAY_MODULES)}")
    for comp in COMPARISONS:
        df = res[(res["comparison.label"] == comp) & (res["sample"] == SAMPLE)]
        if len(df) == 0:
            continue
        print(f"  PLOT: {comp[:50]}...")
        fig = plot_volcano_overlay(df, modules, OVERLAY_MODULES, comp)
        save_figure(fig, str(dir_overlay / f"volcano_cytoskeleton_adhesion_{slugify(comp)}"))

    # ------- 3. Per-module -------
    print("\n=== Per-module ===")
    for comp in COMPARISONS:
        df = res[(res["comparison.label"] == comp) & (res["sample"] == SAMPLE)]
        if len(df) == 0:
            continue

        for i, mod_name in enumerate(all_modules):
            mod_color = COLORS4[i % len(COLORS4)]
            mod_genes = set(modules[modules["module"] == mod_name]["Gene_norm"].unique())

            # Skip if no hits/candidates in this module
            n_hc = ((df["Gene_norm"].isin(mod_genes)) &
                    (df["hit_annotation"].isin(["hit", "candidate"]))).sum()
            if n_hc == 0:
                continue

            fig = plot_volcano_per_module(df, modules, mod_name, comp,
                                          mod_color=mod_color, label_genes=True)
            save_figure(fig, str(dir_per_module / f"volcano_{slugify(mod_name)}_{slugify(comp)}"))

        print(f"  Completed: {comp[:50]}...")

    print(f"\nDone!")
    print(f"  Hits:       {dir_hits}")
    print(f"  Overlays:   {dir_overlay}")
    print(f"  Per-module: {dir_per_module}")


if __name__ == "__main__":
    main()