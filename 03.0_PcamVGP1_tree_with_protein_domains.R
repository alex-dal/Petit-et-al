library(pals)
library(ggplot2)
library(ggtree)
library(treeio)
library(tidyverse)
library(ggplot2)
library(ggnewscale)
library(dplyr)
library(colorspace)
library(tibble)


setwd("/Users/alexandradallaire/genomics/pcam/blast_hits_nt")
tree <- read.tree("PcamV-GP1_Penicillium_hits_1000_top2_nonoverlapping_merged.fa.fasttre")

setwd("/Users/alexandradallaire/genomics/pcam/blast_hits_nt/dna_ori")
tree <- read.tree("PcamV-GP1_Penicillium_hits_6000_top2_nonoverlapping_merged.fa.fasttre") # only 6kb+
tree <- read.tree("O3yPbqOoOPbDYMY-VXyJMA_newick.txt") # only 6kb+ rerooted
tree <- read.tree("boxMXFNZkyyRrWUnDUE-Hw_newick.txt") # only 6kb+ rerooted + no branch length

setwd("/Users/alexandradallaire/genomics/pcam/blast_hits_nt")
domains <- read.delim("PcamV-GP1_Penicillium_hits_1000_top2_nonoverlapping_merged_100longorfs.cdsearchdomains")
hits <- read.delim("PcamV-GP1_Penicillium_hits_description.txt")

# Check tree
p <- ggtree(tree) +
  theme_tree2() 
p
ggsave(filename = "tree.pdf", plot = p, width = 2, height = 11, units = "in", dpi = 600)


# Check tree 6kb
p <- ggtree(tree) +
  geom_tiplab(size=1, align = TRUE, linetype = "blank") +
  xlim(0, 1) +
  theme_tree2() 
p
ggsave(filename = "tree_6kb.pdf", plot = p, width = 9, height = 11, units = "in", dpi = 600)

p <- ggtree(tree) +
  theme_tree2() 
p

# Clean domain and hit names
domains <- domains %>%
  mutate(Query = sub("_orf[0-9]+$", "", Query))
domains <- domains %>%
  mutate(Query = sub("^Q#[0-9]+ - >", "", Query))
hits$Description <- gsub(",", "", hits$Description)
hits$Description <- gsub("\\(", "", hits$Description)
hits$Description <- gsub("\\)", "", hits$Description)
hits$Description <- gsub(":", "_", hits$Description)
hits <- hits %>%
  mutate(Query_Cover = as.numeric(str_remove(Query_Cover, "%")))
colnames(hits)[1] <- "label"

all(domains$Query %in% tree$Query)
setdiff(domains$Query, tree$tip.label)

colnames(domains)[1] <- "label"

domains$label <- factor(domains$label, levels = tree$tip.label)

# Get tree tip positions
tree_data <- p$data
tree_data$label <- gsub("_R_", "", tree_data$label)
# Only tips
tips <- tree_data %>%
  filter(isTip) %>%
  select(label, y)

# Join domain data with tree y coordinates
domains_plot <- left_join(domains, tips, by=c("label"))

# Merge Gal4 and Gal4 supermafily
domains_plot$Short.name <- gsub("GAL4 superfamily", "GAL4", domains_plot$Short.name)
# Remove Tra5 and Glmu superfamily as there is just one and it overlaps with another
domains_plot <- domains_plot[!grepl("Tra5 superfamily", domains_plot$Short.name), ]
domains_plot <- domains_plot[!grepl("glmU superfamily", domains_plot$Short.name), ]
domains_plot <- domains_plot[!grepl("RNase_H_like superfamily", domains_plot$Short.name), ]
domains_plot <- domains_plot[!grepl("rve superfamily", domains_plot$Short.name), ]
domains_plot <- domains_plot[!grepl("Herpes_env superfamily", domains_plot$Short.name), ] # Low p-value, wrong frame


p2 <- p +
  geom_rect(
    data = domains_plot,
    aes(
      xmin = From,
      xmax = To,
      ymin = y - 0.5,
      ymax = y + 0.5,
      fill = Short.name
    ),
    inherit.aes = FALSE     # <- critical
  ) +
  labs(x = "Amino Acid Position", fill = "Domain Type") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  scale_fill_brewer(palette = "Set3")

p2

ggsave(
  filename = "domains.pdf",
  plot = p2,
  width = 4,
  height = 11,
  units = "in"
)

# Bubble plot in one column

# Remove everything before Penicillium in label col
# Example with a data frame called hits or df
hits$label <- sub(".*?(Penicillium.*)", "\\1", hits$label)

# 1. Get numeric y for hits by joining with p$data
tips2 <- tips
tips2$label <- sub(".*?(Penicillium.*)", "\\1", tips2$label)

hits2$Per._ident <- as.numeric(hits2$Per._ident)

hits2 <- hits %>%
  rowwise() %>%
  mutate(
    y = p$data$y[grep(label, p$data$label)][1]  # take first match
  ) %>%
  ungroup()

# 2. Define bubble x position — e.g., 10% beyond max domain end
max_len <- max(domains_plot$To, na.rm = TRUE)
bubble_x <- max_len * 1.1

# 3. Combine plot
p3 <- p +
  geom_rect(
    data = domains_plot,
    aes(
      xmin = From,
      xmax = To,
      ymin = y - 0.5,
      ymax = y + 0.5,
      fill = Short.name
    ),
    inherit.aes = FALSE
  ) +
  geom_point(
    data = hits2,
    aes(
      x = bubble_x,
      y = y,
      size = Query_Cover,
      color = Per._ident
    ),
    inherit.aes = FALSE,
    alpha = 0.7
  ) +
  scale_size_continuous(name = "Query Cover (%)", range = c(0, 1)) +
  scale_color_gradient(low = "blue", high = "red", name = "Percent Identity") +
  labs(x = "Amino Acid Position", fill = "Domain Type") +
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  scale_fill_brewer(palette = "Set3") +
  coord_cartesian(xlim = c(0, bubble_x + max_len*0.1))  # extend x axis to show bubbles

p3

ggsave(
  filename = "domains_bubble.pdf",
  plot = p3,
  width = 5,
  height = 11,
  units = "in"
)


# Bubble plot in two columns
# Assume hits2 already has y coordinates from the tree
max_len <- max(domains_plot$To, na.rm = TRUE)
bubble_x <- max_len * 1.1  # fixed x-position for the circles
offset <- max_len * 0.1  # distance between the two columns
offset2 <- 2 * offset  # position for third column

clade_lookup <- tribble(
  ~Species, ~Section,
  
  # Section Citrina
  "Penicillium_citrinum", "Citrina",
  "Penicillium_manginii", "Citrina",
  "Penicillium_odoratum", "Citrina",
  
  # Section Lanata-Divaricata
  "Penicillium_janthinellum", "Lanata-Divaricata",
  "Penicillium_griseoroseum", "Lanata-Divaricata",
  "Penicillium_brasilianum", "Lanata-Divaricata",
  "Penicillium_daleae", "Lanata-Divaricata",
  
  # Section Aspergilloides
  "Penicillium_antarcticum", "Aspergilloides",
  "Penicillium_frequentans", "Aspergilloides",
  
  # Section Canescentia
  "Penicillium_canescens", "Canescentia",
  "Penicillium_nucicola", "Canescentia",
  
  # Section Brevicompacta
  "Penicillium_brevicompactum", "Brevicompacta",
  "Penicillium_olsonii", "Brevicompacta",
  "Penicillium_salamii", "Brevicompacta",
  
  # Section Ramosa
  "Penicillium_capsulatum", "Ramosa",
  
  # Section Robsamsonia
  "Penicillium_robsamsonii", "Robsamsonia",
  "Penicillium_dipodomyicola", "Robsamsonia",
  
  # Section Roquefortorum
  "Penicillium_roqueforti", "Roquefortorum",
  "Penicillium_paneum", "Roquefortorum",
  
  # Section Osmophila
  "Penicillium_concentricum", "Osmophila",
  "Penicillium_samsonianum", "Osmophila",
  
  # Section Chrysogena
  "Penicillium_chrysogenum", "Chrysogena",
  "Penicillium_rubens", "Chrysogena",
  "Penicillium_nalgiovense", "Chrysogena",

  # Section Fasciculata - camemberti
  "Penicillium_camemberti", "Fasciculata-camemberti",
  "Penicillium_caseifulvum", "Fasciculata-camemberti",
  "Penicillium_commune", "Fasciculata-camemberti",
  "Penicillium_biforme", "Fasciculata-camemberti",
  "Penicillium_palitans", "Fasciculata-camemberti",
  "Penicillium_solitum", "Fasciculata-camemberti",
  "Penicillium_fuscoglaucum", "Fasciculata-camemberti",
  "Penicillium_crustosum", "Fasciculata-camemberti",
  
  # Section Fasciculata - viridicata
  "Penicillium_freii", "Fasciculata-viridicata",
  "Penicillium_polonicum", "Fasciculata-viridicata",

  # Section Penicillium
  "Penicillium_expansum", "Penicillium",
  "Penicillium_italicum", "Penicillium",
  
  # Section Fasciculata - verrucosa
  "Penicillium_verrucosum", "Fasciculata-verrucosa",
  "Penicillium_nordicum", "Fasciculata-verrucosa",
  "Penicillium_thymicola", "Fasciculata-verrucosa",
  
  # Uncertain / unresolved
  "Penicillium_cf.", "Unknown",
  "Penicillium_sp.", "Unknown",
  "Penicillium_mexicanum", "Unknown",
  "TPA_asm:_Penicillium", "Unknown"
)

hits2 <- hits2 %>%
  left_join(clade_lookup, by = "Species")

colnames(hits2)[11] <- "Section"
n_sections <- length(unique(na.omit(hits2$Section)))

# Generate 12 bright, distinct colours
set.seed(123) # Optional
species_cols <- sample(colors(), 15)
species_cols <- as.vector(pals::polychrome(25))[1:10]

#hits2 <- hits2 %>%
#  filter(!is.na(y)) %>%
#  filter(!is.na(Section))

p4 <- p +
  # Domains
  geom_rect(
    data = domains_plot,
    aes(
      xmin = From,
      xmax = To,
      ymin = y - 0.5,
      ymax = y + 0.5,
      fill = Short.name
    ),
    inherit.aes = FALSE
  ) +
  scale_fill_brewer(palette = "Set3") +
  
  # Query_Cover circles
  geom_point(
    data = hits2,
    aes(
      x = bubble_x,
      y = y,
      color = Query_Cover
    ),
    shape = 16,
    size = 1.2,
    inherit.aes = FALSE
  ) +
  scale_color_gradient(low = "gold", high = "blueviolet", name = "Query Cover (%)") +
  
  ggnewscale::new_scale_color() +
  
  # Per._ident circles (second column)
  geom_point(
    data = hits2,
    aes(
      x = bubble_x + offset,
      y = y,
      color = Per._ident
    ),
    shape = 16,
    size = 1.2,
    inherit.aes = FALSE
  ) +
  scale_color_gradient(low = "gold", high = "blueviolet", name = "Percent Identity") +
  
  ggnewscale::new_scale_color() +
  
  # Species column (categorical colors)
  geom_point(
    data = hits2,
    aes(x = bubble_x + offset2, y = y, color = Section),
    size = 1,
    inherit.aes = FALSE
  ) +
  scale_color_manual(
    values = species_cols,
    name = "Penicillium Section"
  )  +
  
  theme_minimal() +
  theme(panel.grid = element_blank()) +
  labs(x = "Amino Acid Position", fill = "Domain Type") +
  coord_cartesian(
    xlim = c(0, bubble_x + offset + offset2 + max_len*0.05)
  )

p4

ggsave(
  filename = "domains_bubbles3_6kb.pdf",
  plot = p4,
  width = 5,
  height = 11,
  units = "in"
)

