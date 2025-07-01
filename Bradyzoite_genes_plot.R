# (Run DEseq as in volcano_plot.R code file)
# Define the genes to plot
genes_to_plot <- c(
  "TGME49_200385",  # BFD1
  "TGME49_311100",  # BFD2
  "TGME49_286090",  # Elf1.2
  "TGME49_251740",  # AP2XI9
  "TGME49_259020",  # BAG1
  "TGME49_291040"   # LDH2
)

# Subset `dge` for the genes to plot
dge_subset <- dge[gene_id %in% genes_to_plot]

# Rename trait factor levels for display
dge_subset[, trait := factor(trait, levels = c(
  "S48_vs_ME49_Tachyzoite",
  "S48_vs_ME49_Brady2",
  "S48_vs_ME49_Brady4",
  "S48_vs_ME49_Brady7"
), labels = c(
  "Tachyzoite",
  "Differentiation Day 2",
  "Differentiation Day 4",
  "Differentiation Day 7"
))]

# Define gene labels
gene_labels <- c(
  "TGME49_200385" = "BFD1",
  "TGME49_311100" = "BFD2",
  "TGME49_286090" = "Elf1.2",
  "TGME49_251740" = "AP2XI9",
  "TGME49_259020" = "BAG1",
  "TGME49_291040" = "LDH2"
)

# Plot

gg <- ggplot(data = dge_subset, aes(x = trait, y = log2FoldChange, fill = trait)) +
  geom_col() +
  geom_hline(yintercept = 0, colour = "grey30") +
  coord_flip() +
  facet_wrap(~gene_id, scales = "free_y", labeller = labeller(gene_id = gene_labels)) +
  theme_light() +  
  scale_fill_manual(values = c(
    "Tachyzoite" = "#0A014F",
    "Differentiation Day 2" = "#ffcad4",
    "Differentiation Day 4" = "#ff97b7",
    "Differentiation Day 7" = "#ff5d8f"
  )) +
  labs(
    title = "Log2 Fold Change of Bradyzoite Markers",
    x = "Contrast",
    y = "Log2 Fold Change"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.text.y = element_text(size = 18),
    axis.title = element_text(size = 18),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "none",  # This line removes the legend
    strip.background = element_rect(fill = "#383F51"),
    strip.text = element_text(color = "white", face = "bold", size = 18)
  )



# Show plot
print(gg)

# Save plot
ggsave("Differentiation_genes.png", gg, width = 20, height = 10, units = "in", dpi = 300)

