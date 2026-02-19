# Phylogenetics and structural analysis

Phylogenetic tree construction pipelines and AlphaFold3 docking configuration for gene family analyses in the paper. These scripts document the computational workflows run on HPC clusters and the AlphaFold3 server; they do not require local data files (input sequences are retrieved from public databases).

## Scripts

| Script | Figure | Description |
|--------|--------|-------------|
| `gene_tree.decarboxylase_pipeline.txt` | Fig 3F | Aromatic amino acid decarboxylase gene tree: BLAST/PROST homolog retrieval → MAFFT alignment → IQ-TREE maximum likelihood tree |
| `gene_tree.transproter_pipeline.txt` | Fig 3G | MFS vesicular transporter gene tree: same pipeline as decarboxylases |
| `gene_tree.gpcr_pipeline.txt` | Fig 5A, S18–S21 | GPCR gene trees for aminergic receptor families |
| `AF3DockingScript.json` | Fig 5A–B | Example AlphaFold3 input configuration for NPY1R-tyramine docking prediction |
| `Halichondria_insilicoPCR_pipeline.txt` | Fig S26 | *Halichondria pinaza* species identification: in silico PCR → 28S rRNA alignment → ML tree |
| `Eunapius_GetOrganelle_pipeline.txt` | Fig S27 | *Eunapius fragilis* species identification: GetOrganelle mitogenome assembly → COI alignment → ML tree |

## Software requirements

These pipelines were run on HPC systems. Key tools (not included in the repository environments):

- **BLAST+** — sequence homolog retrieval
- **PROST** — remote homolog detection via protein structure
- **MAFFT** — multiple sequence alignment
- **IQ-TREE** — maximum likelihood phylogenetic inference
- **AlphaFold3** — protein–ligand docking (accessed via server)
- **GetOrganelle** — organelle genome assembly
- **STAR** — read alignment (referenced in species ID pipelines)

## Notes

- The pipeline `.txt` files document command-line workflows with variable placeholders (e.g., `${SPECIES}`, `${OUTDIR}`). Adapt paths and SLURM parameters for your computing environment.
- Proteome databases referenced in the gene tree pipelines are available from UniProt and NCBI and are too large to include in this repository.
- The `AF3DockingScript.json` is an example input for the AlphaFold3 web server, not an executable script.
