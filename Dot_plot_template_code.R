library(data.table)
library(ggplot2)
library(viridis)
library(limma)
library(DESeq2)
# (Run DEseq as in volcano_plot.R code file)
# Create contrasts
contrasts <- limma::makeContrasts(
 # S48_vs_ME49_Tachyzoite = (S48.Tachyzoite - ME49.Tachyzoite),
  S48_vs_ME49_Brady2 = (S48.Brady_Day_2 - ME49.Brady_Day_2),
  S48_vs_ME49_Brady4 = (S48.Brady_Day_4 - ME49.Brady_Day_4),
  S48_vs_ME49_Brady7 = (S48.Brady_Day_7 - ME49.Brady_Day_7),
  levels = design
)


# Loop to shrink data and merge go terms
merged_list <- vector("list", length = ncol(contrasts))
names(merged_list) <- colnames(contrasts)

for(cntr in names(merged_list)) {
  res_dt <- as.data.table(
    as.data.frame(
      lfcShrink(deseq, type = "ashr", contrast = contrasts[,cntr], quiet = TRUE)
    ),
    keep.rownames = "gene_id"
  )
  res_dt[, trait := cntr]
  
  merged <- merge(res_dt, go_terms, 
                  by.x = "gene_id", by.y = "GID",
                  all.x = TRUE, allow.cartesian = TRUE
  )[!is.na(Name) & Name != "n/a"]
  
  merged_list[[cntr]] <- merged[!duplicated(Name)]
}

merged_data <- rbindlist(merged_list)

# Subset SRS/SAG genes AND filter for >1 Log2fold
dge_subset <- merged_data[
  grepl("^(SRS|SAG)", Name, ignore.case = TRUE) & abs(log2FoldChange) > 1
]

dge_subset[, trait := factor(trait, levels = colnames(contrasts))]

# Rename contrast labels for plotting
dge_subset[, trait := factor(trait, levels = c(
  "S48_vs_ME49_Tachyzoite",
  "S48_vs_ME49_Brady2",
  "S48_vs_ME49_Brady4",
  "S48_vs_ME49_Brady7"
))]

# Define Names
contrast_labels <- c(
  S48_vs_ME49_Tachyzoite = "Tachyzoite",
  S48_vs_ME49_Brady2     = "Differentiation Day 2",
  S48_vs_ME49_Brady4     = "Differentiation Day 4",
  S48_vs_ME49_Brady7     = "Differentiation Day 7"
)


# Plot Dot plot

dot_plot <- ggplot(dge_subset, aes(x = Name, y = log2FoldChange, 
                                   size = abs(log2FoldChange), color = log2FoldChange)) +
  geom_hline(yintercept = 0, color = "black", size = 0.3) + # bolds the intercept 0 line
  geom_point(alpha = 0.8) + # transparency of gene dots
  scale_color_viridis(option = "inferno") +
  scale_size(range = c(3, 7)) +  # dots size range
  facet_grid(rows = vars(trait), scales = "free_y", space = "free_y",
             labeller = labeller(trait = contrast_labels)) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20), # Centered, bold title
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 16),
    axis.title = element_text(size = 18),
    strip.text.y = element_text(size = 20),
    panel.spacing.y = unit(1.2, "lines"), # Spacing between facet panels!
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  ) +
  labs(
    title = "Significant Log2FC of SRS Genes",
    x = "Gene Name",
    y = expression(Log[2]~Fold~Change),
    color = expression(Log[2]~Fold~Change),
    size = expression("Log[2]~Fold~Change")  
  )

ggsave("SRS_Genes_DotPlot_FacetGrid_notachy.pdf", dot_plot, width = 22, height = 20, units = "in")


# Create a new workbook
wb <- createWorkbook()

# Add each contrast's top genes as a separate worksheet
for (cntr in names(merged_list)) {
  addWorksheet(wb, sheetName = cntr)
  writeData(wb, sheet = cntr, merged_list_list[[cntr]])
}

# Save the workbook
saveWorkbook(wb, file = "Top_Genes_by_Contrast.xlsx", overwrite = TRUE)




