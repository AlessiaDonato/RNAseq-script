library(ggplot2)
library(ggrepel)
library(dplyr)

# Load DEG (Differentially Expressed Genes) results table
s1 <- read.delim("res_s1_wt_vs_s1_ko_gene_name.txt")

# Assign colors based on log2FoldChange and adjusted p-value thresholds using dplyr::case_when for readability
# Make sure to change the column names (log2FoldChange, padj) according to your DEG table
s1 <- s1 %>%
  mutate(
    col = case_when(
      abs(log2FoldChange) > 1 & padj < 0.05 ~ "red",           # Significant and large fold change
      abs(log2FoldChange) > 1 & padj >= 0.05 ~ "forestgreen",  # Large fold change but not significant
      abs(log2FoldChange) <= 1 & padj < 0.05 ~ "blue",         # Significant but small fold change
      TRUE ~ "black"                                            # Not meeting any condition (default color)
    )
  )

# Load gene names to annotate on the volcano plot
s1_name <- read.delim("res_s1_wt_vs_s1_ko_gene_name-cg.txt")

# Generate the volcano plot
# - Points colored based on significance and fold change
# - Labels added with geom_text_repel to avoid overlapping text
# Remember to adjust the column names (log2FoldChange, padj, SYMBOL) according to your files
ggplot(s1, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = s1$col) +
  # ylim(c(0, -log10(10e-12))) # Uncomment this line to set y-axis limits if needed
  geom_text_repel(
    data = s1_name,
    aes(x = log2FoldChange, y = -log10(padj), label = SYMBOL),
    max.overlaps = Inf,            # Show all labels without removing any
    box.padding = 0.5,             # Padding around labels to reduce overlap
    point.padding = 0.5,           # Padding around points for label placement
    segment.color = 'grey50'       # Color of lines connecting labels to points
  )
