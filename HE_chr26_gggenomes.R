library(rtracklayer)
library(dplyr)
library(forcats)
library(gggenomes)

gff_files_dir <- "4.2_chr26_HCE_n_flankings_gff3_final/"

gff_files <- paste0(gff_files_dir, c(
#  "Ref_v2_chr26_HCE_and_flankings.gff3",
  "Ref_v3_chr26_HCE_and_flankings.gff3",
  "new_ref_hap1_chr26_HCE_and_flankings.gff3",
#  "new_ref_hap2_chr26_HCE_and_flankings.gff3",
  "BS1_hap1_chr26_HCE_and_flankings.gff3",
  "BS1_hap2_chr26_HCE_and_flankings.gff3",
  "BS2_hap1_chr26_HCE_and_flankings.gff3",
#  "BS2_hap2_chr26_HCE_and_flankings.gff3",
  "BS3_hap1_chr26_HCE_and_flankings.gff3",
  "BS3_hap2_chr26_HCE_and_flankings.gff3",
  "BS4_hap1_chr26_HCE_and_flankings.gff3",
  "BS4_hap2_chr26_HCE_and_flankings.gff3",
  "BS5_hap1_chr26_HCE_and_flankings.gff3",
  "BS5_hap2_chr26_HCE_and_flankings.gff3",
  "BS6_hap1_chr26_HCE_and_flankings.gff3",
  "BS6_hap2_chr26_HCE_and_flankings.gff3",
  "CS2_hap1_chr26_HCE_and_flankings.gff3",
  "CS2_hap2_chr26_HCE_and_flankings.gff3",
  "CS4_hap1_chr26_HCE_and_flankings.gff3",
  "CS4_hap2_chr26_HCE_and_flankings.gff3",
  "CS5_hap1_chr26_HCE_and_flankings.gff3",
  "CS5_hap2_chr26_HCE_and_flankings.gff3",
  "CS7_hap1_chr26_HCE_and_flankings.gff3",
  "CS7_hap2_chr26_HCE_and_flankings.gff3",
  "CS8_hap1_chr26_HCE_and_flankings.gff3",
  "CS8_hap2_chr26_HCE_and_flankings.gff3",
  "CS10_hap1_chr26_HCE_and_flankings.gff3",
  "CS10_hap2_chr26_HCE_and_flankings.gff3",
  "NSSH2_hap1_chr26_HCE_and_flankings.gff3",
  "NSSH2_hap2_chr26_HCE_and_flankings.gff3",
  "NSSH10_hap1_chr26_HCE_and_flankings.gff3",
  "NSSH10_hap2_chr26_HCE_and_flankings.gff3",
  "sprat_chr1_synteny.gff3"
))

# Function to import, filter, modify, and convert to data.frame
process_gff_df <- function(file) {
  gff <- import(file, format = "gff3")
  gff_gene <- gff[gff$type == "gene"]
  gff_gene$type <- "CDS"
  df <- as.data.frame(gff_gene)
  colnames(df)[1] <- "seq_id"
  df$source_file <- basename(file)
  return(df)
}

# Process all files and bind into one data.frame
gff_combined_df <- do.call(rbind, lapply(gff_files, process_gff_df))

gff_combined_df <- gff_combined_df %>%
  mutate(status = ifelse(grepl("null|disrupted", Name, ignore.case = TRUE), "disrupted allele", Target)) %>%
  mutate(seq_id = fct_recode(seq_id,
                             "Ref_v3_hap1" = "Ref_v3_chr26",
                             "Ref_v3_hap2" = "new_ref_hap1_chr26",
                             "Sprat" = "OY735285.1")) %>%
  mutate(seq_id = fct_relabel(seq_id, ~ gsub("_chr26", "", .x))) %>%
  mutate(Target = ifelse(Target == "dnaja3b", "DNAJA3B", ifelse(Target == "hsp70", "HSP70", Target))) %>%
  mutate(status = ifelse(status == "dnaja3b", "DNAJA3B", ifelse(status == "hsp70", "HSP70", status)))


rect_data_Atl <- data.frame(
  ymin = 1.5,
  ymax = 17.5
)

rect_data_Bal <- data.frame(
  ymin = 17.5,
  ymax = 30.5
)

#-------------------------------------------------------------------------------

p <- gggenomes(gff_combined_df) +
  geom_seq() +
  geom_rect(data = rect_data_Atl, aes(ymin = ymin, ymax = ymax), xmin = -50000, xmax = -5000, fill = "#0080FF", alpha = 0.25) +
  geom_rect(data = rect_data_Bal, aes(ymin = ymin, ymax = ymax), xmin = -50000, xmax = -5000, fill = "red", alpha = 0.25) +
  geom_gene(aes(color = Target, fill = status), alpha = 1, size = 2) +
  scale_color_manual(breaks = c("DNAJA3B", "HSP70", "HCE", "CCNB1"), values = c("DNAJA3B" = "#4AC14A", "HSP70" = "#F7BE00", "HCE" = "#0055FF", "CCNB1" = "#FF1212")) +
  scale_fill_manual(breaks = c("DNAJA3B", "HSP70", "HCE", "CCNB1"), ## legend only for these 4
                    values = c("DNAJA3B" = "#7CCD7C", "HSP70" = "#F7C600", "HCE" = "#0066FF", "CCNB1" = "brown1", "disrupted allele" = "white")) + ## "disrupted allele" still used in the plot
  geom_bin_label(size = 4, fontface = "plain") +
  theme(legend.position = "inside",
        legend.position.inside =  c(.75, .25),
        legend.key.size = unit(18, "pt"),
        legend.key.spacing.y = unit(7, "pt"),
        legend.title = element_blank(),
        legend.text = element_text(size = 14, face = "italic"),
        panel.grid = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 14)
  ) +
  guides(color = "none", fill = guide_legend(override.aes = list(color = c("#4AC14A", "#F7BE00", "#0055FF", "#FF1212")), ncol = 1))

pdf("HCE_chr26.pdf", width = 10, height = 7)
print(p)
dev.off()

