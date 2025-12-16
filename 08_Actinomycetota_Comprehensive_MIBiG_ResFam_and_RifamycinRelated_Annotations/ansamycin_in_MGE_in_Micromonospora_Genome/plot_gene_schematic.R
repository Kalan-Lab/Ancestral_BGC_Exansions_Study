#!/usr/bin/env Rscript

library(gggenes)
library(ggplot2)
library(dplyr)
library(readr)

# Configuration
INPUT_FILE <- "genes_for_plotting.tsv"
OUTPUT_FILE <- "gene_schematic.png"

# Color scheme for manual annotations
annotation_colors <- c(
  "AHBA synthesis" = "#FA8072",        # Salmon (light red)
  "PKS" = "#C00000",                   # Red
  "NRPS(-like)" = "#BDC3C7",           # Gray (same as Other)
  "Cytochrome P450" = "#BDC3C7",       # Gray (same as Other)
  "Cytochrome 450" = "#BDC3C7",        # Gray (same as Other)
  "Transporter" = "#3498DB",           # Blue
  "Transposase / Putative Transposase" = "#5DADE2", # Light blue
  "Reverse transcriptase" = "#2E86C1", # Medium blue
  "Reverse transcriptase " = "#2E86C1",# Medium blue (with space)
  "Putative phage integrase" = "#1B4F72", # Dark blue
  "Tailoring enzyme" = "#F39C12",      # Orange/Yellow
  "Regulator" = "#2ECC71",             # Green
  "Other" = "#BDC3C7"                  # Gray
)

# Read the data
cat("Reading gene data...\n")
genes <- read_tsv(INPUT_FILE, show_col_types = FALSE)

# Simplify categories - combine grey categories into "Other" for cleaner legend
genes <- genes %>%
  mutate(
    display_category = case_when(
      manual_annotation %in% c("Cytochrome P450", "Cytochrome 450", "NRPS(-like)") ~ "Other",
      TRUE ~ manual_annotation
    ),
    direction = ifelse(strand == "+", 1, -1)
  )

cat(sprintf("Loaded %d genes\n", nrow(genes)))

# Create the plot
cat("Creating gene schematic with gggenes...\n")

# Calculate midpoints for labels
genes <- genes %>%
  mutate(midpoint = (start + end) / 2)

p <- ggplot(genes, 
            aes(xmin = start, xmax = end, y = molecule, 
                fill = display_category,
                forward = direction)) +
  geom_gene_arrow(arrowhead_height = unit(5, "mm"), 
                  arrowhead_width = unit(2, "mm"),
                  arrow_body_height = unit(5, "mm")) +
  # Add conservation indicators (black stars for conserved genes)
  geom_point(data = filter(genes, conserved == "Yes"),
             aes(x = midpoint, y = 1.12),
             shape = 8, size = 3.5, color = "black", stroke = 1.2,
             inherit.aes = FALSE) +
  scale_fill_manual(values = annotation_colors, 
                    name = "Gene Function") +
  scale_x_continuous(labels = scales::comma) +
  theme_genes() +
  theme(legend.position = "bottom") +
  labs(
    title = sprintf("Biosynthetic Gene Cluster: %s to %s", 
                    genes$locus_tag[1], 
                    genes$locus_tag[nrow(genes)]),
    subtitle = "★ Black star indicates conserved gene",
    x = "Genomic Position (bp)",
    y = NULL
  ) +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE))

# Save the plot
ggsave(OUTPUT_FILE, plot = p, width = 20, height = 6, dpi = 300, bg = "white")
cat(sprintf("\nGene schematic saved to: %s\n", OUTPUT_FILE))

# Print statistics
cat("\nPlot Statistics:\n")
cat(sprintf("  Total genes: %d\n", nrow(genes)))
cat(sprintf("  Forward strand: %d\n", sum(genes$strand == "+")))
cat(sprintf("  Reverse strand: %d\n", sum(genes$strand == "-")))
cat(sprintf("  Conserved genes: %d\n", sum(genes$conserved == "Yes")))
cat(sprintf("  Genomic range: %s - %s bp\n", 
            format(min(genes$start), big.mark = ","),
            format(max(genes$end), big.mark = ",")))

# Print functional categories (using display categories for cleaner output)
cat("\nFunctional Categories (displayed):\n")
func_summary <- genes %>%
  group_by(display_category) %>%
  summarise(count = n(), .groups = 'drop') %>%
  arrange(desc(count))

for (i in 1:nrow(func_summary)) {
  cat(sprintf("  %-35s: %2d gene(s)\n", 
              func_summary$display_category[i], 
              func_summary$count[i]))
}

cat("\tDone!\n")

