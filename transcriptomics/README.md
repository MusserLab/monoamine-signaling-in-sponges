# Transcriptomics

RNA-seq differential expression analysis of tryptamine-treated *Spongilla lacustris* at 1 h, 4 h, and 24 h recovery timepoints. This directory contains the DESeq2 pipeline, GO enrichment analysis, module heatmap generation, and the read mapping documentation.

## Scripts

| Script | Language | Figure | Description |
|--------|----------|--------|-------------|
| `Transcriptomic_analysis_v2.Rmd` | R | Methods | DESeq2 interaction model (`~ Timepoint + Treatment + Timepoint:Treatment`): read count normalization, differential expression, hit/candidate calling |
| `TranscriptomicGOTermAnalysis.Rmd` | R | Fig 4F | GO term enrichment using goseq with EggNOG-based GO annotations |
| `heatmap_transcriptomics.py` | Python | Fig 4G, S17 | Module heatmap of significant transcriptomics genes organized by functional category |
| `transcriptomic_mapping_pipline.txt` | — | Methods | Read mapping pipeline documentation: fastp (QC/trimming) → STAR (alignment) → featureCounts (quantification) |

## Data files

| File | Description |
|------|-------------|
| `data/mapped_transcript_counts.txt` | featureCounts output: read counts per gene across all samples. Input to the DESeq2 pipeline |
| `data/Master_DE_Results_Interaction_V2.tsv` | Pre-computed DESeq2 results for all comparisons (original output without gene short names — names are joined at runtime from `spongilla_gene_names_final.tsv`) |
| `data/transcriptomics_module_assignments.csv` | Significant genes with curated module/category assignments for the heatmap |

| `data/sl_go_eggnog.txt` | Gene-to-GO term mapping (eggNOG mapper + manual corrections; originally `2020_08_13_merged_trinity_gene_gene_ontology_assignments_jake_manual_additions2_forsl2_final2.txt` from Musser et al. 2021) |

Shared reference files used by scripts in this directory (located in `../data/`):
- `spongilla_gene_names_final.tsv` — gene name lookup (join on `Trinity_geneID`, label with `Zang_et_al_2026`)
- `GO_GeneLengths.tsv` — gene lengths for goseq

## Key analytical choices

- **DESeq2 interaction model** tests for treatment effects that differ across timepoints
- **Hit thresholds** (transcriptomics): padj < 0.05 and |log2FC| > 2
- **Candidate thresholds** (transcriptomics): padj < 0.05 and |log2FC| > 0.58
- These thresholds differ from the phosphoproteomics thresholds because of the different dynamic ranges in the two data types

## Running

```r
# DESeq2 pipeline (generates differential expression results)
rmarkdown::render("transcriptomics/Transcriptomic_analysis_v2.Rmd")

# GO enrichment (requires DESeq2 results)
rmarkdown::render("transcriptomics/TranscriptomicGOTermAnalysis.Rmd")
```

```bash
# Heatmap (requires transcriptomics_module_assignments.csv + gene names)
conda activate monoamine-sponges
python transcriptomics/heatmap_transcriptomics.py
```
