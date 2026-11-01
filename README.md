# Petit et al. — analysis scripts

Scripts used to detect, characterise and visualise insertions of the **PcamV-GP1** element across *Penicillium* genomes, annotate the genes and protein domains flanking each insertion, and identify the **Starships** (giant fungal transposons) that carry them.

> **Citation:** Petit *et al.* (under review).

---

## Overview

| Step | Script | What it does |
|---|---|---|
| 1 | `01.0_PcamVGP1_blast.sh` | Builds a combined genome database, finds PcamV-GP1 copies with BLASTn, extracts and aligns the hits, builds an ML tree, and prepares ORFs for domain annotation |
| 2 | `02.0_PcamVGP1_flanks.sh` | Extracts ±50 kb around each insertion and annotates protein domains in the flanks |
| 2a | `02.1_get_orfs.py` | Six-frame translation; writes ORFs ≥ *N* aa with their DNA coordinates in the header |
| 2b | `02.2_clean_cdsearch_outputs.py` | Cleans NCBI CD-Search outputs and parses strand, frame and coordinates from query names |
| 2c | `02.3_merge_cdsearch_domains.py` | Merges overlapping hits of the same domain, keeps the longest domain when different domains overlap, and converts protein positions to genomic coordinates |
| 3 | `03.0_PcamVGP1_tree_with_protein_domains.R` | Plots the PcamV-GP1 tree with domain architectures, BLAST query cover and identity, and *Penicillium* section for each tip |
| 4 | `04.0_PcamVGP1_flanks.R` | Draws gene/domain maps of insertion flanks (Fig. S2) and of Starships (Fig. 6) with **gggenomes** |
| 5 | `05.0_Penicillium_gene_annotation.sh` | Annotates genes in the *Penicillium* genomes with miniprot protein hints + AUGUSTUS (SLURM array job) |
| 6 | `06.0_Starfish_workflow.sh` | Runs the **starfish** workflow to find captain tyrosine recombinases and Starship elements, then groups them into families and haplotypes |
| 7 | `07.0_Starship_orf_annotation.sh` | Extracts Starships that contain PcamV-GP1 and annotates their ORFs and domains (same approach as step 2) |

```
Genomes ──► 01 BLAST PcamV-GP1 ──► hits.bed ──► alignment + tree ──► 03 tree/domain figure
                     │
                     ├──► 02 ±50 kb flanks ──► ORFs ──► CD-Search ──► 04 flank maps (Fig. S2)
                     │
Genomes ──► 05 gene annotation ──► 06 starfish ──► Starships ──► 07 ORFs/CD-Search ──► 04 Starship maps (Fig. 6)
```

---

## Requirements

**Command-line tools**

- [BLAST+](https://blast.ncbi.nlm.nih.gov/) (`makeblastdb`, `blastn`)
- [bedtools](https://bedtools.readthedocs.io/)
- [samtools](http://www.htslib.org/)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [bioawk](https://github.com/lh3/bioawk)
- [MAFFT](https://mafft.cbrc.jp/alignment/software/)
- [IQ-TREE 3](http://www.iqtree.org/) (`iqtree3`)
- [miniprot](https://github.com/lh3/miniprot)
- [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus)
- [starfish](https://github.com/egluckthaler/starfish), plus its dependencies (HMMER, MMseqs2, MCL, MUMmer…)

**Python ≥ 3.8**

- `biopython`
- `pandas`

**R ≥ 4.2**

- `ggtree`, `treeio` (Bioconductor)
- `ggplot2`, `tidyverse`, `dplyr`, `tibble`
- `ggnewscale`
- `gggenomes`
- `colorspace`
- `pals`

**Web service**

- [NCBI Batch CD-Search](https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi), default settings

---

## Inputs

- Genome assemblies (`.fna`), one directory per genome. NCBI assemblies go in `GCA_*` / `GCF_*` directories. Other assemblies have their directory name added to each FASTA header so contig names stay unique.
- `PcamV-GP1.fa`: the PcamV-GP1 query sequence.
- `protein_evidence.faa`: protein evidence for miniprot/AUGUSTUS.
- `genomes_augustus.txt`: soft-masked genome paths for the SLURM array (`*.rm.fa`), one per line.
- `PcamVGP1_starships.bed`: coordinates of Starships carrying PcamV-GP1, selected from the starfish output.

---

## Usage

The scripts are **step-by-step analysis notebooks**, not a turnkey pipeline. Run them section by section from a working directory that holds the inputs above.

### 1. Find PcamV-GP1 copies and build the tree

```bash
bash 01.0_PcamVGP1_blast.sh
```

Key settings:

- BLASTn e-value ≤ 1e-5, hits > 2 kb.
- Six-frame ORFs ≥ 100 aa.
- IQ-TREE with 1000 UFBoot and 1000 SH-aLRT replicates.

Manual steps:

- **Extend the early-diverging hits.** In `hits.bed`, extend these hits by 6,434 bp (strand-aware) to recover truncated domains. Save the edited file as `hits.bed2`.
- **Add the query sequence.** Paste the PcamV-GP1 query sequence into `PcamV-GP1_hits_extracted.fna` before alignment.
- **Run CD-Search.** Submit `PcamV-GP1_hits_extracted_100longorfs.faa` to NCBI Batch CD-Search.

### 2. Annotate the insertion flanks

```bash
bash 02.0_PcamVGP1_flanks.sh
```

This extracts ±50 kb around each hit and translates ORFs with `02.1_get_orfs.py`. It then splits the ORFs into batches of 999 sequences for CD-Search; submit each batch and save the results as `cdsearch_batches/*.txt`. Finally, it cleans and merges the domain hits and writes `cdsearch_outs_merged.bed` along with `PcamVGP1_loci_flanks_fixhead.fasta`.

### 3–4. Figures

Run in R or RStudio:

- `03.0_PcamVGP1_tree_with_protein_domains.R`: tree + domains + query cover, identity and section (`domains_bubbles3_6kb.pdf`)
- `04.0_PcamVGP1_flanks.R`: `Starship_Figure_S2.pdf` (flanks) and `Starship_Figure_6.pdf` (Starships)

### 5. Gene annotation (SLURM)

```bash
sbatch --array=1-$(wc -l < genomes_augustus.txt) 05.0_Penicillium_gene_annotation.sh
```

This uses miniprot protein alignments as `CDSpart` hints for AUGUSTUS (`--species=anidulans`, soft-masking on).

### 6. Starship detection with starfish

`06.0_Starfish_workflow.sh` follows the [starfish tutorial](https://github.com/egluckthaler/starfish/wiki). It runs these steps:

- **Gene finder:** `annotate` → `consolidate` → `sketch`
- **Element finder:** `insert` (two rounds, the second with `--pid 80 --hsp 750`) → `summarize` → YR superfamily assignment → `sim`/`group` into families and haplotypes → `pair-viz`
- **Flanks:** `flank`

Several steps take 1–2 h each, and `pair-viz` takes much longer.

### 7. Annotate Starships carrying PcamV-GP1

```bash
bash 07.0_Starship_orf_annotation.sh
```

This works the same way as step 2, but runs on `PcamVGP1_starships.fasta` and writes `starship_cdsearch_outs_merged_filtered.bed` and `PcamVGP1_starships_fixhead.fasta`.

---

## Notes and known quirks

- **Helper script names.** The `.sh` scripts call the helpers without their numeric prefix (`get_orfs.py`, `clean_cdsearch_outputs.py`, `merge_cdsearch_domains.py`). Rename or symlink them, or edit the calls.
- **`02.1_get_orfs.py` settings are hardcoded.** The input/output files are currently set for the Starship run (step 7) with `min_aa = 100`. For the flank run (step 2), switch to the commented `PcamVGP1_loci_flanks*` names and set `min_aa = 300`.
- **`02.3_merge_cdsearch_domains.py` output names.** It processes every `*_cln.txt` in the current directory and writes `*_merged_filtered.txt`. Step 2 expects `cdsearch_outs_merged.txt`, so rename one or the other.
- **The R scripts were written for interactive use.** Update the `setwd()` paths and pick the `read.tree()` line for the tree you want. Some objects, such as `hits2`, are defined further down than where they're first used, so run the blocks in order.
- **`06.0_Starfish_workflow.sh` can't be run directly.** It's a lab-notebook record of the commands and contains non-shell lines (`//`, `----`). Run its commands one at a time.
- **`01.0` uses `cat */*.fna`.** For non-NCBI genomes this picks up both the original and the `_prefixed.fna` files. Move or remove the originals first, or narrow the glob.
- **CD-Search is a manual web step.** Remove the `#` comment header lines from the downloaded results before cleaning them.

---

## Repository structure

```
.
├── 01.0_PcamVGP1_blast.sh
├── 02.0_PcamVGP1_flanks.sh
├── 02.1_get_orfs.py
├── 02.2_clean_cdsearch_outputs.py
├── 02.3_merge_cdsearch_domains.py
├── 03.0_PcamVGP1_tree_with_protein_domains.R
├── 04.0_PcamVGP1_flanks.R
├── 05.0_Penicillium_gene_annotation.sh
├── 06.0_Starfish_workflow.sh
├── 07.0_Starship_orf_annotation.sh
└── README.md
```

## Contact

Alexandra Dallaire alexandra.dallaire@bcm.ulaval.ca — open an issue in this repository for questions.

