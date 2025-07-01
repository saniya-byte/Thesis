# Load packages
library(tximport)
library(rhdf5)
library(tidyverse)
library(biomaRt)
library(pheatmap)
library(plotly)
library(clusterProfiler)
library(DESeq2)
library(RColorBrewer)
library(org.Bt.eg.db)
library(org.Mm.eg.db)
library(ggrepel)
library(dplyr)
library(DESeq2)
library(data.table)
library(limma)
library(ashr)
library(ggplot2)
library(ggrepel)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnhancedVolcano)
library(openxlsx)  # for Excel output


# Set working directory
setwd("/Volumes/Saniya's Space/RNA_seq/Kallisto_output_ME49 copy")

# Load sample table
sample_table <- read.csv("run.table.csv")

# Convert sample_table to data.table
sample_table <- as.data.table(sample_table)

# Subset samples with condition not equal to 'HFF'
#sample_table <- sample_table[sample_table$condition != 'HFF', ]
# Replace all white space with underscore in the condition column
sample_table[, condition := gsub(' ', '_', condition)]

# Create file paths for homan and specific transcriptomes
homan_files <- paste0(
  pull(sample_table, "sample_name"),
  "/kallisto_output_homan/abundance.h5"
)

# Human gene mapping
homan_gene_map <- read.csv("Human_transcript_gene_mapping.txt", 
                           header = FALSE, 
                           sep = "\t", 
                           col.names = c("Transcript", "GeneID"))

# Import kallisto output for both transcriptomes
homan_data <- tximport(files = homan_files,
                       type = "kallisto",
                       tx2gene = homan_gene_map,
                       ignoreAfterBar = TRUE)

# Get GO terms - create go_terms
# Extract gene IDs from the homan_gene_map
gene_ids <- unique(homan_gene_map$GeneID)
# Remove version numbers
gene_ids_clean <- sub("\\..*", "", gene_ids)

# Query the database for GO terms and their descriptions
go_terms <- AnnotationDbi::select(org.Hs.eg.db, 
                                  keys = gene_ids_clean, 
                                  columns = c("GO", "SYMBOL"), 
                                  keytype = "ENSEMBL")
# Convert to data.table
go_terms <- as.data.table(go_terms)

go_terms <- go_terms[, .(
  GO = paste(GO, collapse = ','),
  Name = SYMBOL[1]
), by = ENSEMBL]


# Add a new 'group' column in sample_table for DESeq design
sample_table[, group := paste(strain, condition, sep = '.')]
design <- model.matrix(~0 + group, sample_table)
colnames(design) <- sub('^group', '', colnames(design))

# DESeq analysis
deseq <- DESeqDataSetFromTximport(homan_data, colData = sample_table, design = design)
deseq <- DESeq(deseq)

# Variance stabilising transformation
vsd <- vst(assay(deseq))
colnames(vsd) <- sample_table$sample_name

# Select the 500 most variable genes
select <- order(apply(vsd, 1, var), decreasing = TRUE)[1:500]

# PCA calculation
pca <- prcomp(t(vsd[select, ]))
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
d <- data.table(PC1 = pca$x[, 1], PC2 = pca$x[, 2], sample_name = colnames(vsd))
pcaout <- merge(d, sample_table, by = 'sample_name')

# Rename condition labels for PCA plot
pcaout$condition <- factor(pcaout$condition, 
                           levels = c("Tachyzoite", "Brady_Day_2", "Brady_Day_4", "Brady_Day_7"),
                           labels = c("Tachyzoite", "Bradyzoite Day 2", "Bradyzoite Day 4", "Bradyzoite Day 7"))

# Define time point colors
time_colors <- c("Tachyzoite" = "#ED254E", 
                 "Bradyzoite Day 2" = "#F9DC5C", 
                 "Bradyzoite Day 4" = "#C2EABD", 
                 "Bradyzoite Day 7" = "#011936")

# Define strain shapes
strain_shapes <- c("RH" = 17,   # Triangle
                   "S48" = 16,  # Circle
                   "ME49" = 15) # Square

# PCA Plot 
gg <- ggplot(data = pcaout, aes(x = PC1, y = PC2, 
                                colour = condition, shape = strain)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = time_colors) +  
  scale_shape_manual(values = strain_shapes) +  
  xlab(sprintf('PC1 (%.1f%% variance explained)', percentVar[1] * 100)) +
  ylab(sprintf('PC2 (%.1f%% variance explained)', percentVar[2] * 100)) +
  ggtitle('Principal Component Analysis of Samples') +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_line(color = "gray90", linetype = "dashed"),
        legend.position = "right")

print(gg)
# ggsave("Human_PCA_final.pdf", width = 8, height = 6)

# Define contrasts for differential expression analysis
contrasts <- limma::makeContrasts(
  S48_vs_ME49_Tachyzoite = (S48.Tachyzoite) - (ME49.Tachyzoite),
  S48tachy_RHtachy       = (S48.Tachyzoite) - (RH.Tachyzoite),
  S48_specific           = (S48.Tachyzoite) - ((ME49.Tachyzoite) + (RH.Tachyzoite))/2,
  # bradyzoite Day 2 comparisons
  Tachyzoite_to_Bradyzoite = (S48.Tachyzoite - S48.Brady_Day_2) - (ME49.Tachyzoite - ME49.Brady_Day_2),
  S48_vs_ME49_Diff_Day2    = S48.Brady_Day_2 - ME49.Brady_Day_2,
  # bradyzoite Day 4 comparisons
  S48_vs_ME49_Diff_Day4    = S48.Brady_Day_4 - ME49.Brady_Day_4,
  Day2_to_Day4             = (S48.Brady_Day_2 - S48.Brady_Day_4) - (ME49.Brady_Day_2 - ME49.Brady_Day_4),
  # bradyzoite Day 7 comparisons
  S48_vs_ME49_Diff_Day7    = S48.Brady_Day_7 - ME49.Brady_Day_7,
  Day4_to_Day7             = (S48.Brady_Day_4 - S48.Brady_Day_7) - (ME49.Brady_Day_4 - ME49.Brady_Day_7),
  levels = design
)

# Ensure that the contrasts matrix and design matrix match
if (!identical(rownames(contrasts), colnames(design))) {
  stop("Mismatch between contrast row names and design column names!")
}

# Initialise list for differential expression results
dge <- list()

# Loop for differential expression contrasts
for (cntr in colnames(contrasts)) {
  message("Processing contrast: ", cntr)
  tryCatch({
    res <- lfcShrink(deseq, type = 'ashr', contrast = contrasts[, cntr], quiet = TRUE)
    res <- as.data.table(as.data.frame(res), keep.rownames = 'ENSEMBL')
    res[, trait := cntr]
    dge[[length(dge) + 1]] <- res
  }, error = function(e) {
    message("Error processing contrast ", cntr, ": ", e$message)
  })
}

# Combine all results into a single data table
dge <- rbindlist(dge, fill = TRUE)
# Remove version numbers in ENSEMBL IDs
dge$ENSEMBL <- sub("\\..*", "", dge$ENSEMBL)

# Merge differential expression results with GO terms
merged_data <- merge(dge, go_terms, by = "ENSEMBL", all.x = TRUE, allow.cartesian = TRUE)


### Volcano Plot per Contrast with p-value cutoff 0.001

# Loop through each contrast and generate volcano plots
for (contrast_name in valid_contrasts) {
  message("Generating volcano plot for contrast: ", contrast_name)
  
  # Subset data for the current contrast
  contrast_data <- merged_data[trait == contrast_name & !is.na(Name) & Name != "n/a"]
  
  if (nrow(contrast_data) == 0) {
    message("No data available for contrast: ", contrast_name)
    next
  }
  
  # Remove duplicate gene names
  contrast_data_unique <- contrast_data[!duplicated(Name), ]
  
  # Define consistent thresholds (same as EnhancedVolcano)
  p_cutoff <- 0.0001
  fc_cutoff <- 3.0
  
  # Select genes to label based on the exact volcano thresholds
  genes_to_label <- contrast_data_unique[
    pvalue < p_cutoff & abs(log2FoldChange) > fc_cutoff, Name]
  
  # Generate plot
  volcano_plot <- EnhancedVolcano(
    contrast_data_unique,
    lab = contrast_data_unique$Name,
    selectLab = genes_to_label,
    x = 'log2FoldChange',
    y = 'pvalue',
    xlab = bquote(~Log[2]~ 'fold change'),
    pCutoff = p_cutoff,         
    FCcutoff = fc_cutoff,
    pointSize = 3.0,
    labSize = 4.0,
    labCol = 'black',
    labFace = 'bold',
    boxedLabels = FALSE,
    colAlpha = 0.5,
    legendPosition = 'right',
    legendLabSize = 12,
    legendIconSize = 3.0,
    drawConnectors = TRUE,
    widthConnectors = 0.5,
    colConnectors = 'grey50',
    title = contrast_name,
    col = c("#FFD166", "#06D6A0", "#EF476F", "#26547C"),
    axisLabSize = 12,
    gridlines.major = FALSE,
    gridlines.minor = FALSE,
    max.overlaps = 20  # Limit the number of label overlaps to declutter
  )
  
  # Save the volcano plot
  pdf_filename <- paste0("Human_vol_", contrast_name, ".pdf")
  ggsave(pdf_filename, plot = volcano_plot, device = "pdf", width = 10, height = 8)
  
  print(volcano_plot)
}


### Save DEseq contrast data to Excel workbook 

# Create a new Excel workbook

wb <- createWorkbook()

# Loop through each unique contrast to add filtered data as a separate sheet
for (contrast in valid_contrasts) {
  # Subset data for the current contrast
  contrast_data <- merged_data[trait == contrast]
  
  # Filter for genes with log2FoldChange >= 1 and pvalue < 0.001
  contrast_data_filtered <- contrast_data[abs(log2FoldChange) >= 1 & pvalue < 0.001, ]
  

  sheet_name <- substr(contrast, 1, 31)
  addWorksheet(wb, sheetName = sheet_name)
  
  # Write the filtered contrast data into the sheet
  writeData(wb, sheet = sheet_name, x = contrast_data_filtered)
}

# Save the workbook
saveWorkbook(wb, file = "All_Contrasts.xlsx", overwrite = TRUE)



### Extract and Save top genes per contrast as a CSV file

extract_top_genes <- function(merged_data) {
  # Filter significant genes using p-value cutoff 0.001 and fold change > 2
  significant_genes <- merged_data[pvalue < 0.001 & abs(log2FoldChange) > 2, ]
  
  # Select top 10 upregulated and downregulated genes per contrast
  top_upregulated <- significant_genes[order(-log2FoldChange), .SD[1:10], by = trait]
  top_downregulated <- significant_genes[order(log2FoldChange), .SD[1:10], by = trait]
  
  # Combine results
  top_genes <- rbind(top_upregulated, top_downregulated)
  
  # Keep only relevant columns
  top_genes <- top_genes[, .(Name, trait, log2FoldChange, ENSEMBL)]
  
  # Save results to CSV
  fwrite(top_genes, "Human_top_genes_per_contrast.csv")
  
  print(top_genes)
}

extract_top_genes(merged_data)



