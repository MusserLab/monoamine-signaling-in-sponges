# Microscopy

Quantitative image analysis scripts for canal kymograph measurements and immunofluorescence quantification.

## Scripts

| Script | Language | Figure | Description |
|--------|----------|--------|-------------|
| `tyramine_sl_kymograph_analysis_final.Rmd` | R | Fig 1K, M | Kymograph analysis: reads canal diameter measurements over time, computes Pearson distance between paired canals to quantify coordination of contraction/relaxation across DMSO and tyramine conditions |
| `vinc_analysis_count_batch.ijm` | ImageJ macro | Fig 5L | Batch counting of vinculin foci in immunofluorescence images |
| `vinc_analysis_areaFoci_batch.ijm` | ImageJ macro | Fig 5M | Batch measurement of vinculin foci area |
| `Pmyl_intensity_analysis_batch_v2.ijm` | ImageJ macro | Fig 5N | Batch measurement of phospho-myosin light chain (pMYL) fluorescence intensity on F-actin structures |

## Data files

| File | Description |
|------|-------------|
| `data/kymo_metadata.tsv` | Kymograph metadata: sample IDs, conditions (DMSO/tyramine), replicate numbers |
| `data/DMSO_R*_kym.csv` | Kymograph measurements for DMSO-treated replicates (6 files) |
| `data/Tyr_R*_kym.csv` | Kymograph measurements for tyramine-treated replicates (7 files) |

## Running

```r
# Kymograph analysis
rmarkdown::render("microscopy/tyramine_sl_kymograph_analysis_final.Rmd")
```

The ImageJ macros (`.ijm` files) are run within [FIJI/ImageJ](https://fiji.sc/). They use dialog-based file/folder selection — no hardcoded paths. Open the macro in FIJI via Plugins → Macros → Run, then select your input image directory when prompted.
