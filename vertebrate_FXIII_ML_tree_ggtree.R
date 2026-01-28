library(tidyverse)
library(ggtree)
library(ggtreeExtra)
library(treeio)
library(ape)
library(cowplot)
library(ggnewscale)

FTG_tree <- read.newick("AA_based_ML_tree_rerooted_ape_edited.nwk", node.label='support')
FTG_tree@phylo$tip.label

## Rename some tip labels
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Clupea_harengus_FTG_chr17_1"] <- "Clupea_harengus_FTG1 (chr17)"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Clupea_harengus_FTG_chr17_2"] <- "Clupea_harengus_FTG2 (chr17)"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Clupea_harengus_FTG_chr17_3"] <- "Clupea_harengus_FTG3 (chr17)"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Clupea_harengus_FTG_chr15_1"] <- "Clupea_harengus_FTG4 (chr15)"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Clupea_harengus_FTG_chr15_2"] <- "Clupea_harengus_FTG5 (chr15)"

FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Clupea_harengus_f13a1b_chr18"] <- "Clupea_harengus_f13a1b (chr18)"

FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Danio_rerio_FTG_1"] <- "Danio_rerio_FTG1"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Danio_rerio_FTG_2"] <- "Danio_rerio_FTG2"

FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Oncorhynchus_mykiss_FTG_1"] <- "Oncorhynchus_mykiss_FTG1"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Oncorhynchus_mykiss_FTG_2"] <- "Oncorhynchus_mykiss_FTG2"
FTG_tree@phylo$tip.label[FTG_tree@phylo$tip.label == "Oncorhynchus_mykiss_FTG_3"] <- "Oncorhynchus_mykiss_FTG3"


metadata <- data.frame(tip_labels = FTG_tree@phylo$tip.label) %>%
  mutate(category = ifelse(grepl("Clupea_harengus", tip_labels), "Clupea_harengus", "other"))

clades_df <- data.frame(clade = c("Teleost FTG", "Teleost F13A"), node = c(50, 72))


p1 <- ggtree(FTG_tree, linetype = NA) %<+% metadata +
  geom_highlight(
    data = subset(clades_df, clade == "Teleost FTG"),
    aes(node = node, fill = clade),
    extend = 0.345, alpha = 1, type = "gradient", gradient.direction = "tr"
  ) +
  geom_highlight(
    data = subset(clades_df, clade == "Teleost F13A"),
    aes(node = node, fill = clade),
    extend = 0.385, alpha = 1, type = "gradient", gradient.direction = "tr"
  ) +

  # —— Branch coloring by branch length (continuous) ——
  geom_tree(aes(color = branch.length), linewidth = 0.7) +
  scale_color_viridis_c(option = "plasma", name = "Substitutions per site") +
  
  # Start a NEW color scale for tip labels so it doesn’t clash with branch colors
  ggnewscale::new_scale_color() +
  
  geom_tiplab(aes(color = category), fontface = "plain", align = FALSE, size = 3.5,
              linetype = "dotted", linesize = 0.3, hjust = 0, offset = 0) +
  
  scale_color_manual(values = c("Clupea_harengus" = "#CC0000", "other" = "black"), guide = "none") +
  
  scale_fill_manual(values = c("Teleost FTG" = "#FFD5FF", "Teleost F13A" = "#D1D1FF"), breaks = c("Teleost FTG", "Teleost F13A"), name = NULL) +
  geom_treescale(y = -0.7, linesize = 0.6, width = NULL, offset = 0.5, fontsize = 5) +
  xlim(0, 1.9) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(.84, .33),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 14, face = "plain"),
    legend.key.size = unit(20, "pt"),
    plot.margin = margin(5, 0, 5, 0, "pt")
  ) +
  
  guides(
    colour = guide_colorbar(order = 1),  # branch length first
    fill   = guide_legend(order = 2)     # clade second
  )
p2 <- flip(p1, 15, 16)
p3 <- flip(p2, 13, 14)
p4 <- flip(p3, 79, 82)

ggsave("ML_tree_rooted_on_arowana.pdf", plot = p4, device = "pdf", width = 10, height = 10)

