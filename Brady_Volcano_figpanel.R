# Load necessary libraries
library(EnhancedVolcano)
library(patchwork)
library(data.table)
library(DESeq2)
library(limma)
# ... (other required libraries)

# (Run DEseq and create design model matric as in volcano_plot.R code file)

# create contrasts
contrasts <- limma::makeContrasts(
  # bradyzoite Day 2 comparisons
  Tachyzoite_to_Bradyzoite = (S48.Tachyzoite - S48.Brady_Day_2) - (ME49.Tachyzoite - ME49.Brady_Day_2),
  S48_vs_ME49_Diff_Day2    = S48.Brady_Day_2 - ME49.Brady_Day_2,
  #bradyzoite Day 4 comparisons
  S48_vs_ME49_Diff_Day4    = S48.Brady_Day_4 - ME49.Brady_Day_4,
  Day2_to_Day4             = (S48.Brady_Day_2 - S48.Brady_Day_4) - (ME49.Brady_Day_2 - ME49.Brady_Day_4),
  #bradyzoite Day 7 comparisons
  S48_vs_ME49_Diff_Day7    = S48.Brady_Day_7 - ME49.Brady_Day_7,
  Day4_to_Day7             = (S48.Brady_Day_4 - S48.Brady_Day_7) - (ME49.Brady_Day_4 - ME49.Brady_Day_7),
  levels = design
)

# Create list to store results and volcano plots
dge_list <- list()
selected_genes_list <- list()
volcano_plots <- list()

# Loop over contrasts
for(cntr in colnames(contrasts)) {
  
  message("Processing contrast: ", cntr)
  
  # Shrink log2 fold changes for the current contrast using ashr
  res <- lfcShrink(deseq, type = 'ashr', contrast = contrasts[, cntr], quiet = TRUE)
  
  # Convert results to data.table and keep gene IDs as a column
  res <- as.data.table(as.data.frame(res), keep.rownames = 'gene_id')
  
  # Add a column indicating the contrast
  res[, trait := cntr]
  
  # Save full result
  dge_list[[cntr]] <- res
  
  # Merge with GO term annotations 
  merged_res <- merge(res, go_terms, by.x = "gene_id", by.y = "GID", all.x = TRUE, allow.cartesian = TRUE)
  
  # Remove rows with missing gene names and remove duplicates
  merged_res <- merged_res[!is.na(Name) & Name != "n/a", ]
  merged_unique <- merged_res[!duplicated(Name), ]
  
  # Select top 10 up-regulated and top 10 down-regulated genes
  top_genes <- merged_unique[order(log2FoldChange, decreasing = TRUE)][1:10, ]
  bottom_genes <- merged_unique[order(log2FoldChange)][1:10, ]
  selected_genes <- rbind(top_genes, bottom_genes)
  selected_genes[, Contrast := cntr]
  selected_genes_list[[cntr]] <- selected_genes
  
  # Create label vector: label only the selected genes
  top10_up <- merged_unique[order(log2FoldChange, decreasing = TRUE)][1:10, Name]
  top10_down <- merged_unique[order(log2FoldChange)][1:10, Name]
  labels <- ifelse(merged_unique$Name %in% c(top10_up, top10_down), merged_unique$Name, NA)
  
  # Generate the volcano plot
  volcano_plot <- EnhancedVolcano(
    merged_unique,
    lab = labels,
    x = 'log2FoldChange',
    y = 'pvalue',
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = 10e-3,  #  p-value 0.001
    FCcutoff = 1, #  log2FC filtering
    pointSize = 4.0,
    labSize = 4.5,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = FALSE,
    colAlpha = 0.5,
    legendPosition = 'right',
    legendLabSize = 14,
    legendIconSize = 4.0,
    drawConnectors = TRUE,
    widthConnectors = 1.0,
    colConnectors = 'black',
    title = cntr,
    col = c("#ffbe0b", "#fb5607", "#ff006e", "#8338ec"),
    axisLabSize = 14,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    max.overlaps = Inf
  )
  
  # Print the plot to the active device
  print(volcano_plot)
  
  # Store the plot in a list
  volcano_plots[[cntr]] <- volcano_plot
  
  # save to file if needed
}


### Create bradyzoite volacno panel figure


library(patchwork)
library(ggplot2)
library(grid)

# Reduce margins for each volcano plot to make them all fit on one page.
# Adjust the unit values as needed.
volcano_plots <- lapply(volcano_plots, function(p) {
  p + theme(
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "lines")
  )
})


# Combine plots by bradyzoite stage comparisons
brady_day2_panel <- volcano_plots[["Tachyzoite_to_Bradyzoite"]] + 
  volcano_plots[["S48_vs_ME49_Diff_Day2"]]
brady_day4_panel <- volcano_plots[["S48_vs_ME49_Diff_Day4"]] + 
  volcano_plots[["Day2_to_Day4"]]
brady_day7_panel <- volcano_plots[["S48_vs_ME49_Diff_Day7"]] + 
  volcano_plots[["Day4_to_Day7"]]

# Assemble the full panel: increased relative space for plots and minimal white space.
full_panel <- (brady_day2_panel / brady_day4_panel / brady_day7_panel) +
  plot_layout(guides = "collect", heights = c(5, 5, 5)) &
  theme(
    legend.position = "right",
    legend.text = element_text(size = 8),       # smaller legend text
    legend.title = element_text(size = 10),       # smaller legend title
    legend.key.size = unit(0.5, "cm"),            # smaller legend keys
    plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
    plot.spacing = unit(0.2, "cm")                # minimal spacing between panels
  )

# Save as panel fig
ggsave("Bradyzoite_Full_Panel_p001.pdf", full_panel,
       width = 12, height = 16, units = "in", dpi = 600)
