library(dplyr)
library(data.table)
library(reshape2)
library(tidyr)
library(stringr)
library(readxl)
library(Hmisc)
library(readr)
library(dplyr)
library(tidyverse)
library(gggenomes)
library(ggplot2)
library(lang)
library(colorspace)

################################################################################
# Plot ORF annotations surrounding PcamV-GP1 insertions
################################################################################

s1 <- read_seqs("PcamVGP1_loci_flanks_fixhead.fasta")
g1 <- read_feats("cdsearch_outs_merged.bed")

gggenomes(seqs = s1, genes = g1) +
  geom_seq() + 
  geom_gene() +
  geom_gene_text(aes(label = name), angle = 45, size = 0.8, nudge_y = 0.2, hjust = -0.05) +
  geom_bin_label(size=2)

# Colour palette for maximum divergence
all_ogs <- sort(unique(c(
  g1$name)))
og_colours <- setNames(
  qualitative_hcl(
    n = length(all_ogs),
    h = c(0, 422),   # give range, not vector
    c = 100,
    l = 65
  ),
  all_ogs
)  

p <- gggenomes(seqs = s1, genes = g1) +
  geom_seq() + 
  geom_gene(aes(fill = name), colour = NA) +
  scale_fill_manual(values = og_colours, drop = FALSE) +
  geom_gene_text(aes(label = name), angle = 45, size = 1, nudge_y = 0.1, hjust = -0.05) +
  geom_bin_label(size=2) +
  theme(legend.position = "none")

ggsave(
  filename = "Starship_Figure_S2.pdf",
  plot = p,
  width = 10,
  height = 17,
  units = "in"
)


################################################################################
# Plot ORF annotations on Starships annotated by Starfish
################################################################################

s1 <- read_seqs("PcamVGP1_starships_fixhead.fasta")
g1 <- read_feats("starship_cdsearch_outs_merged_filtered.bed")

gggenomes(seqs = s1, genes = g1) +
  geom_seq() + 
  geom_gene() +
  geom_gene_text(aes(label = name), angle = 45, size = 0.8, nudge_y = 0.2, hjust = -0.05) +
  geom_bin_label(size=2)

# Colour palette for maximum divergence
all_ogs <- sort(unique(c(
  g1$name)))
og_colours <- setNames(
  qualitative_hcl(
    n = length(all_ogs),
    h = c(0, 422),   # give range, not vector
    c = 100,
    l = 65
  ),
  all_ogs
)  

p <- gggenomes(seqs = s1, genes = g1) +
  geom_seq() + 
  geom_gene(aes(fill = name), colour = NA) +
  scale_fill_manual(values = og_colours, drop = FALSE) +
  geom_gene_text(aes(label = name), angle = 45, size = 1, nudge_y = 0.1, hjust = -0.05) +
  geom_bin_label(size=2) +
  theme(legend.position = "none")

ggsave(
  filename = "Starship_Figure_6.pdf",
  plot = p,
  width = 30,
  height = 17,
  units = "in"
)
