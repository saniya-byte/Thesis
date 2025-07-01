# Load packages
library(AnnotationDbi)
library(DESeq2)
library(EnhancedVolcano)
library(RColorBrewer)
library(ashr)
library(biomaRt)
library(clusterProfiler)
library(data.table)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(limma)
library(org.Bt.eg.db)
library(org.Mm.eg.db)
library(org.TgondiiME49.eg.db)
library(pheatmap)
library(plotly)
library(rhdf5)
library(tidyverse)
library(tximport)


# Set working directory
setwd("/Volumes/Saniya\'s\ Space/RNA_seq/Kallisto_output_ME49\ copy")

# Load sample table 
sample_table <- read.csv("run.table.csv")

library(data.table)

# Convert sample_table to data.table
sample_table <- as.data.table(sample_table)

#subset with condition ! = not equal
sample_table <- sample_table[sample_table$condition != 'HFF', ]
#replace all white space with underscore
sample_table[, condition := gsub(' ', '_', condition)]

specific_files <- paste0(
  pull(sample_table, "sample_name"),
  "/kallisto_output_specific/abundance.h5"
)

specific_gene_map <- read.csv("ME49_transcript_gene_mapping.txt", header = FALSE, sep = "\t", col.names = c("Transcript", "GeneID"))

specific_data <- tximport(files = specific_files,
                          type = "kallisto",
                          tx2gene = specific_gene_map)

### Gene onotology terms
# Extract gene IDs from the specific_gene_map
gene_ids <- unique(specific_gene_map$GeneID)

# Get GO terms, descriptions and Name from ME49 annotation db using AnnotationDbi package
go_terms <- AnnotationDbi::select(org.TgondiiME49.eg.db, 
                                  keys = gene_ids, 
                                  columns = c("GO", "description", "Name"), 
                                  keytype = "GID")

# Group creates a tag which identifies each condition = combining strain and condition into a single label
sample_table[, group :=paste(strain,condition, sep='.')]

# Construct a design matrix for the DESeq2 model
# ~0 removes the intercept so each group gets its own column 
design <- model.matrix(~0 + group, sample_table)
# remove group prefix from column names
colnames(design) <- sub('^group', '', colnames(design))

# Create a DESeq2 dataset object, by importing transcript quantification data (from tximport)
deseq <- DESeqDataSetFromTximport(specific_data, colData=sample_table, design=design)
deseq <- DESeq(deseq)

# variance-stabilizing transformation applies a transformation to stabilise the variance across genes 
#(which can vary due to differences in sequencing depth or gene expression level).

vsd <- vst(deseq@assays@data$counts)

#Assigns human-readable sample names to the columns of the transformed data (vsd).

colnames(vsd) <- sample_table$sample_name

# Create contrasts to compare how gene expression changes between conditions differ across strains
# Using a full model like this combines log fold changes with statistical testing (p-values),
# while using information from the entire dataset to improve variance estimation and reduce noise
contrasts <- limma::makeContrasts(
  Tachyzoite_to_Bradyzoite = (S48.Tachyzoite - S48.Brady_Day_2)  - (ME49.Tachyzoite - ME49.Brady_Day_2),
  brady2_Brady4 = (S48.Brady_Day_2 - S48.Brady_Day_4)  - (ME49.Brady_Day_2 - ME49.Brady_Day_4),
  levels = design
)


# Check if row names of contrasts match column names of design matrix
stopifnot(identical(rownames(contrasts), colnames(design)))


# Loop over each defined contrast 
# Use DESeq2's lfcShrink to get shrinkage-estimated LFCs for each contrast
# This improves effect size estimates, especially for low-count or noisy genes
dge <- list()
for(cntr in colnames(contrasts)) {
  print(cntr)
  res <- lfcShrink(deseq, type='ashr', contrast=contrasts[, cntr], quiet=TRUE)
  res <- as.data.table(as.data.frame(res), keep.rownames='gene_id')
  res[, trait :=cntr]
  dge[[length(dge) + 1]] <- res
}
dge <- rbindlist(dge)


# Some genes are annotated with multiple GO terms — Collapse GO terms for each gene into one row
go_terms <- as.data.table(go_terms)
go_terms <- go_terms[, .(
  GO = paste(GO, collapse = ','),
  description = description[1],  
  Name = Name[1]  
), by = GID]

# Merge the differential expression results with the GO annotations
# One gene may have multiple GO terms → allow.cartesian = TRUE keeps all combinations 

merged_data <- merge(dge, go_terms, by.x = "gene_id", by.y = "GID", all.x = TRUE, allow.cartesian = TRUE)



# Initialize an empty list for storing results
dge <- list()

### save the contrast data as a tsv file
# Loop through each contrast to process and save results
for(cntr in colnames(contrasts)) {
  print(cntr)
  
  # Perform the shrinkage estimation for the current contrast
  res <- lfcShrink(deseq, type='ashr', contrast=contrasts[, cntr], quiet=TRUE)
  res <- as.data.table(as.data.frame(res), keep.rownames='gene_id')
  res[, trait := cntr]
  
  # Append the result to the list
  dge[[length(dge) + 1]] <- res
}

# Combine all results into a single data table
dge <- rbindlist(dge)

# Merge the differential gene expression results with GO terms
merged_data <- merge(dge, go_terms, by.x = "gene_id", by.y = "GID", all.x = TRUE, allow.cartesian = TRUE)

# Save the merged results as a TSV file for each contrast
for (contrast_name in unique(merged_data$trait)) {
  contrast_data <- merged_data[trait == contrast_name]
  
  # Save the result to a TSV file
  fwrite(contrast_data, paste0("contrast_", contrast_name, "_results.tsv"), sep = "\t")
  

}

### create a volcano plot for each contrast
# Enhanced volano plot loop
for (contrast_name in unique(merged_data$trait)) {
  # Filter data for the specific contrast
  contrast_data <- merged_data[trait == contrast_name]
  
  # Remove rows where 'Name' is 'n/a' and duplicates
  contrast_data <- contrast_data[!is.na(Name) & Name != "n/a", ]
  contrast_data_unique <- contrast_data[!duplicated(Name), ]
  
  # Sort and select top and bottom expressed genes
  top_genes <- contrast_data_unique[order(log2FoldChange, decreasing = TRUE)][1:10, ]
  bottom_genes <- contrast_data_unique[order(log2FoldChange)][1:10, ]
  highlight_genes <- c('BAG1', 'LDH2', 'ENO1', 'ATG-9', 'SRS9', 'ELF1.2', 'AP2IX9', 'BFD1', 'BFD2')
  highlight_unique <- contrast_data_unique[contrast_data_unique$Name %in% highlight_genes, ]
  selected_genes <- rbind(top_genes, bottom_genes, highlight_unique)
  
  # Create a label vector
  labels <- ifelse(contrast_data_unique$Name %in% selected_genes$Name, contrast_data_unique$Name, NA)
  
  # Create the Enhanced Volcano plot
  volcano_plot <- EnhancedVolcano(contrast_data_unique,
                                  lab = labels,
                                  x = 'log2FoldChange',
                                  y = 'pvalue',
                                  xlab = bquote(~Log[2]~ 'fold change'),
                                  pCutoff = 10e-14,
                                  FCcutoff = 2.0,
                                  pointSize = 4.0,
                                  labSize = 4.5,
                                  labCol = 'black',
                                  labFace = 'bold',
                                  boxedLabels = FALSE,
                                  colAlpha = 0.5, # makes the points transparent
                                  legendPosition = 'right',
                                  legendLabSize = 14,
                                  legendIconSize = 4.0,
                                  drawConnectors = TRUE,
                                  widthConnectors = 1.0,
                                  colConnectors = 'black',
                                  title = contrast_name,  # Set the title
                                  col = c("#ffbe0b", "#fb5607", "#ff006e", "#8338ec"),
                                  axisLabSize = 14,
                                  gridlines.major = FALSE,
                                  gridlines.minor = FALSE,
                                  max.overlaps = 50)  
  
  
  # Save the plot as a PDF
  #pdf(paste0("volcano_plot_", contrast_name, ".pdf"), width = 8, height = 6)
  print(volcano_plot)
  volcano_plot
  dev.off()
}

