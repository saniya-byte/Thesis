library(patchwork)

# Define contrasts for the panel figure
contrasts <- limma::makeContrasts(
  S48_vs_ME49_Tachyzoite = S48.Tachyzoite - ME49.Tachyzoite,
  S48tachy_RHtachy       = S48.Tachyzoite - RH.Tachyzoite,
  levels = design
)

# 2) Create a list to merge data
merged_list <- vector("list", length = ncol(contrasts))
names(merged_list) <- colnames(contrasts)

#loop over contrasts to shrink log2fold and merge go terms
for(cntr in names(merged_list)) {
  res_dt <- as.data.table(
    as.data.frame(
      lfcShrink(deseq, type = "ashr", contrast = contrasts[,cntr], quiet = TRUE)
    ),
    keep.rownames = "gene_id"
  )
  res_dt[, trait := cntr]
  
  # GO annotations
  
  merged <- merge(res_dt, go_terms, 
                  by.x = "gene_id", by.y = "GID",
                  all.x = TRUE, allow.cartesian = TRUE
  )[!is.na(Name) & Name != "n/a"]
  
  #remove duplicated gene names
  merged_list[[cntr]] <- merged[!duplicated(Name)]
}

# merge contrasts into one table

merged_data <- rbindlist(merged_list)

# Helper function to build a dot‑plot for any gene family
make_dot_plot <- function(pattern, title_text, fc_threshold = NULL) {
  subset_dt <- merged_data[grepl(pattern, Name, ignore.case = TRUE)]
  
  # Optional, filter by log2fc
  if (!is.null(fc_threshold)) {
    subset_dt <- subset_dt[abs(log2FoldChange) > fc_threshold]
  }
  
  # trait order on y axis
  
  subset_dt[, trait := factor(trait, levels = names(contrast_labels))]
  
  # make ggplot
  
  ggplot(subset_dt, aes(x = Name, y = log2FoldChange,
                        size = abs(log2FoldChange), color = log2FoldChange)) +
    geom_hline(yintercept = 0, color = "black", size = 0.3) +
    geom_point(alpha = 0.8) +
    scale_color_viridis(option = "inferno") +
    scale_size(range = c(3, 7)) +
    facet_grid(rows = vars(trait), scales = "free_y", space = "free_y",
               labeller = labeller(trait = contrast_labels)) +
    theme_minimal(base_size = 16) +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y  = element_text(size = 16),
      axis.title   = element_text(size = 18),
      strip.text.y = element_text(size = 20),
      panel.spacing.y = unit(1.2, "lines"),
      legend.title   = element_text(size = 16),
      legend.text    = element_text(size = 14)
    ) +
    labs(
      title = title_text,
      x     = "Gene Name",
      y     = expression(Log[2]~Fold~Change),
      color = expression(Log[2]~Fold~Change),
      size  = expression("|Log[2]~Fold~Change|")
    )
}

# Build each plot — only SRS/SAG uses fc_threshold = 1
plot_MIC <- make_dot_plot("^MIC",    "Log2FC of MIC Genes")
plot_SRS <- make_dot_plot("^(SRS|SAG)", "Log2FC of SRS/SAG Genes", fc_threshold = 1)
plot_AP2 <- make_dot_plot("^AP2",    "Log2FC of AP2 Genes")

# Stack them vertically
panel <- plot_MIC / plot_SRS / plot_AP2 + plot_layout(ncol = 1, guides = "collect") # collect keys - so only on key per panem

ggsave("Panel_MIC_SRS_AP2_DotPlots_filtered.pdf", panel, width = 18, height = 20, units = "in")
