# An ancient monoaminergic signaling system coordinates contractility in a nerveless sponge

**Rong Xuan Zang<sup>1</sup>\*, Nawaphat Malaiwong<sup>1</sup>&dagger;, Ling Wang<sup>2</sup>&dagger;&Dagger;, Jamie D. Maziarz<sup>1</sup>&dagger;, Kejue Jia<sup>1</sup>&dagger;, Bernhard Drotleff<sup>3</sup>, Frank Stein<sup>4</sup>, Sean Liu<sup>1</sup>, Mah Noor<sup>1</sup>, C. Jackson Roberts<sup>1</sup>, Mandy Rettel<sup>4</sup>, Jennifer J. Schwarz<sup>4</sup>, Aissam Ikmi<sup>5</sup>, Shigeki Watanabe<sup>6,7,8</sup>, Robert Prevedel<sup>2,5</sup>, Noriko Funayama<sup>9</sup>, Michael P. O'Donnell<sup>1,10</sup>, Jacob M. Musser<sup>1,10,11</sup>\***

<sup>1</sup> Department of Molecular, Cellular and Developmental Biology, Yale University; New Haven, 06511, USA
<br><sup>2</sup> Cell Biology and Biophysics Unit, European Molecular Biology Laboratory; Heidelberg, 69117, Germany
<br><sup>3</sup> Metabolomics Core Facility, European Molecular Biology Laboratory; Heidelberg, 69117, Germany
<br><sup>4</sup> Proteomics Core Facility, European Molecular Biology Laboratory; Heidelberg, 69117, Germany
<br><sup>5</sup> Developmental Biology Unit, European Molecular Biology Laboratory; Heidelberg, 69117, Germany
<br><sup>6</sup> Department of Cell Biology, Johns Hopkins University School of Medicine; Baltimore, 21205, USA
<br><sup>7</sup> Solomon H. Snyder Department of Neuroscience, Johns Hopkins University School of Medicine; Baltimore, 21205, USA
<br><sup>8</sup> The Center for Cell Dynamics, Johns Hopkins University School of Medicine; Baltimore, 21205, USA
<br><sup>9</sup> Department of Biophysics, Kyoto University; Kyoto, 606-8502, Japan
<br><sup>10</sup> Wu Tsai Institute, Yale University; New Haven, 06511, USA
<br><sup>11</sup> Systems Biology Institute, Yale West Campus, West Haven, 06516, USA

&dagger; These authors contributed equally to this work.

&Dagger; Current address: College of Optical and Electronic Technology, China Jiliang University, Hangzhou, Zhejiang, 310018, China

\* Corresponding authors: jacob.musser@yale.edu, zangrongxuan@gmail.com

**Abstract:** Chemical neurotransmission is central to coordinating animal physiology and behavior, yet when these signaling systems originated is unclear. Sponges are early-diverging animals that lack neurons and muscles, yet still coordinate contraction and relaxation of their filter-feeding water canals. Here, we show that *Spongilla lacustris* synthesizes the monoamines tryptamine, phenethylamine, and tyramine to regulate distinct canal behaviors. We identify previously uncharacterized decarboxylases and vesicular transporters coexpressed in secretory neuroid and metabolic cell types, and show that tryptamine activates GPCR signaling and Rho GTPases to remodel actomyosin and adhesion complexes in contractile canal epithelia, driving localized constrictions and whole-body deflations. These findings define a monoaminergic signaling system that couples secretory and contractile cells, revealing an ancestral mode of chemical coordination later elaborated in synaptic nervous systems.

---

## Repository overview

This repository contains the data processing pipelines, statistical analyses, and figure-generating code for Zang, Malaiwong, Wang, Maziarz, Jia et al. (2026). Scripts are organized by analysis type and map directly to the paper's main and supplementary figures.

## Figure-to-script mapping

| Figure | Panel | Script | Directory |
|--------|-------|--------|-----------|
| Fig 1 | K, M | `tyramine_sl_kymograph_analysis_final.Rmd` | `microscopy/` |
| Fig 3 | F | `gene_tree.decarboxylase_pipeline.txt` | `phylogenetics_and_structure/` |
| Fig 3 | G | `gene_tree.transproter_pipeline.txt` | `phylogenetics_and_structure/` |
| Fig 4 | Methods | `tryptamine_phosphoproteomics_analysis_cleaned.qmd` | `phosphoproteomics/` |
| Fig 4 | B | `volcano_plot.py` | `phosphoproteomics/` |
| Fig 4 | C | `PhosphoGOTermAnalysis.Rmd` | `phosphoproteomics/` |
| Fig 4 | D | `heatmap_figure_4D.py` | `phosphoproteomics/` |
| Fig 4 | E | `phospho_NO_tryptamine_comparison.qmd`, `NO_tryptamine_scatter.py` | `phosphoproteomics/` |
| Fig 4 | Methods | `Transcriptomic_analysis_v2.Rmd` | `transcriptomics/` |
| Fig 4 | F | `TranscriptomicGOTermAnalysis.Rmd` | `transcriptomics/` |
| Fig 4 | G | `heatmap_transcriptomics.py` | `transcriptomics/` |
| Fig 5 | A | `gene_tree.gpcr_pipeline.txt` | `phylogenetics_and_structure/` |
| Fig 5 | A–B | `AF3DockingScript.json` | `phylogenetics_and_structure/` |
| Fig 5 | L | `vinc_analysis_count_batch.ijm` | `microscopy/` |
| Fig 5 | M | `vinc_analysis_areaFoci_batch.ijm` | `microscopy/` |
| Fig 5 | N | `Pmyl_intensity_analysis_batch_v2.ijm` | `microscopy/` |
| Fig S16 | — | `heatmap_supplemental.py` | `phosphoproteomics/` |
| Fig S17 | — | `heatmap_transcriptomics.py` | `transcriptomics/` |
| Fig S19–S20 | — | `gene_tree.gpcr_pipeline.txt` | `phylogenetics_and_structure/` |
| Fig S19, S21 | — | `AF3DockingScript.json` | `phylogenetics_and_structure/` |
| Fig S27 | — | `Halichondria_insilicoPCR_pipeline.txt` | `phylogenetics_and_structure/` |
| Fig S28 | — | `Eunapius_GetOrganelle_pipeline.txt` | `phylogenetics_and_structure/` |
| Methods | — | `transcriptomic_mapping_pipline.txt` | `transcriptomics/` |

## Directory structure

```
monoamine-signaling-in-sponges/
├── data/                              # Shared reference files (gene names, GO annotations)
├── phosphoproteomics/                 # TMT phosphoproteomics pipeline and figures (7 scripts)
│   └── data/                          #   Co-located phosphoproteomics data files
├── transcriptomics/                   # RNA-seq pipeline and figures (4 scripts)
│   └── data/                          #   Co-located transcriptomics data files
├── phylogenetics_and_structure/       # Phylogenetic trees and AlphaFold3 docking (6 scripts)
├── microscopy/                        # Kymograph analysis and ImageJ macros (4 scripts)
│   └── data/                          #   Kymograph measurement files
├── renv.lock                          # R package versions
├── environment.yml                    # Python/conda environment
└── outs/                              # Generated outputs (not tracked in git)
```

Each directory has its own `README.md` with detailed descriptions of its scripts, data files, and expected outputs.

## Setup

### Prerequisites

- **R** >= 4.3 — [download](https://cran.r-project.org/)
- **renv** (R package manager) — install from an R console with `install.packages("renv")`. renv records the exact package versions used in this project so results are reproducible; see the [renv documentation](https://rstudio.github.io/renv/) for details.
- **conda** (Python environment manager) — install [Miniconda](https://docs.anaconda.com/miniconda/) or [Miniforge](https://github.com/conda-forge/miniforge). conda creates an isolated Python environment with the correct package versions; see the [conda documentation](https://docs.conda.io/) for details.
- **Quarto** >= 1.3 (for rendering `.qmd` scripts) — [install](https://quarto.org/docs/get-started/)
- **Git LFS** (for cloning large data files) — [install](https://git-lfs.com/)

### 1. Clone the repository (with Git LFS)

This repository uses [Git LFS](https://git-lfs.com/) to store large data files (>50 MB). You must install Git LFS before cloning:

```bash
# Install Git LFS (if not already installed)
# macOS:
brew install git-lfs
# Ubuntu/Debian:
sudo apt-get install git-lfs

# Initialize Git LFS
git lfs install

# Clone the repository (LFS files are downloaded automatically)
git clone https://github.com/musserlab/monoamine-signaling-in-sponges.git
cd monoamine-signaling-in-sponges
```

The following files are stored in Git LFS:
- `phosphoproteomics/data/limma_results_annotated.rds` (~91 MB)
- `phosphoproteomics/data/limma_results_annotated.tsv` (~367 MB)
- `phosphoproteomics/data/mdata.rds` (~98 MB)
- `transcriptomics/data/Master_DE_Results_Interaction_V2.tsv` (~71 MB)

If you cloned without Git LFS, run `git lfs pull` to download the large files.

### 2. Set up the R environment

Open R/RStudio/Positron from the repository root directory, then:

```r
# Install all R packages (CRAN + Bioconductor, exact versions pinned in renv.lock)
renv::restore()
```

`renv::restore()` installs all R packages at the exact versions used in the paper, including Bioconductor packages (limma, DESeq2, vsn, etc.).

**Troubleshooting:** If `renv::restore()` fails on specific packages, you can install them individually with `install.packages()` (CRAN) or `BiocManager::install()` (Bioconductor), then retry `renv::restore()` for the rest.

### 3. Set up the Python environment

```bash
conda env create -f environment.yml
conda activate monoamine-sponges
```

This creates a conda environment with numpy, pandas, scipy, matplotlib, seaborn, pyreadr, and adjustText.

### 4. Run scripts

All R scripts use `here::here()` for file paths. Run them from the repository root directory:

```r
# Example: render the core phosphoproteomics pipeline
quarto::quarto_render("phosphoproteomics/tryptamine_phosphoproteomics_analysis_cleaned.qmd")
```

Python scripts use `Path(__file__).resolve().parent` for relative paths and can be run from any directory:

```bash
conda activate monoamine-sponges
python phosphoproteomics/volcano_plot.py
```

## Data availability

| Data type | Repository | Accession |
|-----------|-----------|-----------|
| Raw metabolomics data | Metabolomics Workbench | [PR002929](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Project&ProjectID=PR002929), [ST004638](https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST004638) |
| Raw mass spectrometry | PRIDE | [PXD073078](https://www.ebi.ac.uk/pride/archive/projects/PXD073078) |
| Raw and processed DNA and RNA sequencing datasets | NCBI | [PRJNA1425054](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1425054) |

Processed intermediate data files needed to run the scripts in this repository are included directly (see each directory's `README.md` for descriptions).

## Gene naming

Gene names in this repository follow the conventions described in Zang et al. (2026). The file `data/spongilla_gene_names_final.tsv` maps Trinity transcript IDs to standardized gene names (column `Zang_et_al_2026`), including paralog suffixes (A, B, C) disambiguated by single-cell RNA-seq expression level.

## Citation

If you use code or data from this repository, please cite:

> Zang RX, Malaiwong N, Wang L, Maziarz JD, Jia K, Drotleff B, Stein F, Liu S, Noor M, Roberts CJ, Rettel M, Schwarz JJ, Ikmi A, Watanabe S, Prevedel R, Funayama N, O'Donnell MP, Musser JM. An ancient monoaminergic signaling system coordinates contractility in a nerveless sponge. *bioRxiv* (2026). doi: [10.64898/2026.02.14.704838v1](https://doi.org/10.64898/2026.02.14.704838v1)

## Contributors

- **Rong Xuan Zang** ([@ViciousPlants](https://github.com/ViciousPlants)) — phosphoproteomics and transcriptomics analysis scripts, GO enrichment, heatmaps, volcano plots, phylogenetic pipelines, microscopy analysis
- **Kejue Jia** ([@jkjium](https://github.com/jkjium)) — phylogenetic pipelines and PROST ortholog annotation
- **Jacob M. Musser** ([@JMusser499](https://github.com/JMusser499)) — phosphoproteomics pipeline, NO comparison analysis, repository assembly and documentation

## AI assistance

We utilized ChatGPT 5.2, Gemini 3, and Claude Opus 4.5 and 4.6 to assist with code generation and managing git. All AI-generated code was carefully reviewed and validated by the authors.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.
