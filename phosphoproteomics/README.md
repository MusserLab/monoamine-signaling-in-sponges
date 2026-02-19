# Phosphoproteomics

Quantitative TMT phosphoproteomics analysis of tryptamine-treated *Spongilla lacustris*. This directory contains the core statistical pipeline (FragPipe output to limma differential analysis), GO enrichment, comparison with NO-responsive phosphosites, and all figure-generating scripts for the phosphoproteomics panels.

## Scripts

| Script | Language | Figure | Description |
|--------|----------|--------|-------------|
| `tryptamine_phosphoproteomics_analysis_cleaned.qmd` | R (Quarto) | Fig 4, Methods | Core pipeline: read FragPipe TSVs, batch correction (ComBat), VSN normalization, kNN imputation, control-ratio normalization, limma differential analysis, hit/candidate calling |
| `PhosphoGOTermAnalysis.Rmd` | R | Fig 4C | GO term enrichment using goseq with PROST-based GO annotations. Tests enrichment for each of 5 key comparisons (3 timepoints + 2 difference-in-differences) |
| `phospho_NO_tryptamine_comparison.qmd` | R (Quarto) | Fig 4E, text | Compares NO-responsive phosphosites (2 min treatment) with tryptamine time course. Quantifies overlap, directional concordance, and generates summary statistics |
| `volcano_plot.py` | Python | Fig 4B | Volcano plots showing hits and candidates per comparison, with module-colored overlays and per-module panels |
| `heatmap_figure_4D.py` | Python | Fig 4D | Heatmap of 88 curated genes organized by signaling module, with per-timepoint DMSO-normalized log2 fold changes |
| `heatmap_supplemental.py` | Python | Fig S16 | Supplemental per-module heatmaps for all modules, with column splitting for large modules |
| `NO_tryptamine_scatter.py` | Python | Fig 4E | Scatter plot comparing tryptamine 15 min vs NO 2 min log2 fold changes with Pearson correlation |

## Data files

### Raw inputs (from FragPipe)

| File | Description |
|------|-------------|
| `data/phosphosites_tryptamine_phosphoproteomics.tsv` | FragPipe phosphosite-level quantification (TMT reporter ions) |
| `data/input_tryptamine_phosphoproteomics.tsv` | FragPipe protein-level quantification (input/lysate channel) |
| `data/metadata_tryptamine_phosphoproteomics.csv` | Sample metadata: file paths, TMT channels, conditions, replicates, batches |
| `data/geneID_proteinID_lookup.tsv` | Protein accession to Trinity gene ID mapping |
| `data/spongilla_proteome.fasta` | Predicted proteome (for peptide-to-protein mapping in NO comparison) |

### Pipeline outputs (from `tryptamine_phosphoproteomics_analysis_cleaned.qmd`)

| File | Description |
|------|-------------|
| `data/limma_results_annotated.rds` | Limma differential analysis results with gene annotations (R binary, 91 MB) |
| `data/limma_results_annotated.tsv` | Same as above in tab-separated format for Python scripts (367 MB) |
| `data/mdata.rds` | Tidy long-format measurement data with control-normalized ratios (R binary, 98 MB) |
| `data/conditions.csv` | Sample condition labels |

### Module assignments

| File | Description |
|------|-------------|
| `data/modules_long_RAW_PHOSPHO_RZ_merged.tsv` | Module assignments for raw phospho hits/candidates, with `Plot_Module` and `Label_Name` columns |
| `data/figure_4D_primary_module_assignments.tsv` | Curated 88-gene module assignments for the Fig 4D heatmap |
| `data/gene_list_with_modules_Fig4D_Fig4E.tsv` | Combined module list used by the NO scatter plot |

### NO comparison data

| File | Description |
|------|-------------|
| `data/NO_phosphoproteomics.csv` | Phosphoproteomics from NO-treated *Spongilla* (external dataset) |
| `data/NO_phosphosites_mapped.csv` | NO phosphosites mapped to tryptamine dataset coordinates |
| `data/NO_scatter_label_candidates_RZ.tsv` | Curated gene label candidates for scatter plot variants |

### GO enrichment

| File | Description |
|------|-------------|
| `data/GO_enrichment_background_genes_gene_level.tsv` | Background gene universe for goseq |

Shared reference files used by scripts in this directory (located in `../data/`):
- `spongilla_gene_names_final.tsv` — gene name lookup (join on `Trinity_geneID`, label with `Zang_et_al_2026`)
- `sl_go_merged_id.tsv` — gene-to-GO term mapping
- `GO_GeneLengths.tsv` — gene lengths for goseq length-bias correction

## Key analytical choices

- **Raw phospho** (`sample == "phospho"`) is used for all figures, not input-normalized phospho
- **Per-timepoint DMSO normalization**: all heatmaps and visualizations subtract matched DMSO control at each timepoint, not a single shared DMSO baseline
- **Hit thresholds**: FDR < 0.05 and |log2FC| >= 1
- **Candidate thresholds**: FDR < 0.2 and |log2FC| >= 0.585
- **5 key comparisons**: Tryp vs DMSO at 3, 15, 30 min; difference-in-differences for 3-to-15 and 15-to-30 min transitions

## Running

```r
# Core pipeline (generates limma_results_annotated.rds, mdata.rds, conditions.csv)
quarto::quarto_render("phosphoproteomics/tryptamine_phosphoproteomics_analysis_cleaned.qmd")

# GO enrichment (requires limma results)
rmarkdown::render("phosphoproteomics/PhosphoGOTermAnalysis.Rmd")

# NO comparison (requires limma results + NO data)
quarto::quarto_render("phosphoproteomics/phospho_NO_tryptamine_comparison.qmd")
```

```bash
# Python figure scripts (require limma_results_annotated.tsv + mdata.rds)
conda activate monoamine-sponges
python phosphoproteomics/volcano_plot.py
python phosphoproteomics/heatmap_figure_4D.py
python phosphoproteomics/heatmap_supplemental.py
python phosphoproteomics/NO_tryptamine_scatter.py
```
