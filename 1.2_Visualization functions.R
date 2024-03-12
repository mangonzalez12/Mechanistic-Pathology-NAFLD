#VOLCANO PLOTS
# Assuming list_df is a list of data frames
# Load necessary libraries
library(ggplot2)
library(scales)
library(gridExtra)
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(ComplexUpset)
library(readxl)
library(fgsea)
library(tidyverse)
library(patchwork)
library(ggplotify)
library(umap)
library(grid)
library(plotly)
library(colorspace)
library(RColorBrewer)
library(VennDiagram)
library(DESeq2)
library(ggplot2)
library(utils)
library(readr)
library(dplyr)
library(pheatmap)
library(biomaRt)
library(orca)
library(cowplot)
library(caret)
#library(kaleido)

# Volcano plot
create_volcano_plot <- function(df, plot_title) {
  # Assuming your data frame has columns named "log2FoldChange," "pvalue," "padj," and "Symbol"
  
  # Your volcano plot code here
  p <- ggplot(df, aes(x = log2FoldChange, y = -log10(pvalue))) +
    geom_point(aes(color = case_when(
      padj < 0.01 & log2FoldChange < 0 ~ "#273871",
      padj < 0.01 & log2FoldChange > 0 ~ "#ea4c3b",
      TRUE ~ "grey"
    )), alpha = 0.7) +
    geom_text(
      aes(label = Symbol),
      check_overlap = TRUE, vjust = 1.5, hjust = 1.5, size = 3, color = "black", na.rm = TRUE
    ) +  # Add gene names
    scale_color_identity() +
    labs(
      title = plot_title,
      x = "Log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "blue") +
    annotate("text", x = max(df$log2FoldChange), y = -log10(0.01),
             label = "P-value = 0.01", vjust = -0.5, hjust = 1, color = "blue") +
    theme_bw() +
    # Format y-axis labels without decimals
    scale_y_continuous(labels = scales::number_format(accuracy = 1))
  
  # Return the plot
  return(p)
}



#GSEA function
GSEA = function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  myGO = fgsea::gmtPathways(GO_file)
  
  fgRes <- fgsea::fgsea(pathways = myGO,
                        stats = gene_list,
                        minSize=15, ## minimum gene set size
                        maxSize=400, ## maximum gene set size
                        nperm=10000) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval) %>% 
    arrange(desc(NES))
  message(paste("Number of signficant gene sets =", nrow(fgRes)))
  
  message("Collapsing Pathways -----")
  concise_pathways = collapsePathways(data.table::as.data.table(fgRes),
                                      pathways = myGO,
                                      stats = gene_list)
  fgRes = fgRes[fgRes$pathway %in% concise_pathways$mainPathways, ]
  message(paste("Number of gene sets after collapsing =", nrow(fgRes)))
  
  fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
  filtRes = rbind(head(fgRes, n = 10),
                  tail(fgRes, n = 10 ))
  
  total_up = sum(fgRes$Enrichment == "Up-regulated")
  total_down = sum(fgRes$Enrichment == "Down-regulated")
  header = paste0("Top 10 (Total pathways: Up=", total_up,", Down=",    total_down, ")")
  
  colos = setNames(c("firebrick2", "dodgerblue2"),
                   c("Up-regulated", "Down-regulated"))
  
  g1= ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
    geom_point( aes(fill = Enrichment, size = size), shape=21) +
    scale_fill_manual(values = colos ) +
    scale_size_continuous(range = c(2,10)) +
    geom_hline(yintercept = 0) +
    coord_flip() +
    labs(x="Pathway", y="Normalized Enrichment Score",
         title=header)
  
  output = list("Results" = fgRes, "Plot" = g1)
  return(output)
}

# Plot and save the results plot from GSEA
plot_GSEA <- function(res, title_plot, filename){
  # Make it into a ggplot
  p <- grid.draw(res$Plot)
  # Save file
  ggsave(filename, p,  width = 12, height = 8, dpi = 400)
  print(p)
}



# Run Cluster Profiler
performCluster_profiler <- function(gene_list, pvalue_cutoff = 0.01, qvalue_cutoff = 0.2) {
  require("org.Hs.eg.db")
  require("clusterProfiler")
  
  # Change ENSEMBL id to ENTREZ
  entrez_gene_list <- mapIds(org.Hs.eg.db, keys = gene_list, keytype = "ENSEMBL", column = "ENTREZID")
  
  # Perform gene enrichment analysis using clusterProfiler
  enrich_result <- enrichKEGG(
    gene          = entrez_gene_list,
    organism      = "hsa",
    # keyType       = "ENTREZID",  # Set keyType to "ENSEMBL"
    pvalueCutoff  = pvalue_cutoff,
    qvalueCutoff  = qvalue_cutoff
  )
  # Convert gene ids to gene symbols to make plots more readable
  enrich_result2 <- setReadable(enrich_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  
  return(enrich_result2)
}

### Perform multiple lists
perform_multipleCluster_profiler <- function(gene_lists, pvalue_cutoff = 0.01, qvalue_cutoff = 0.2) {
  require("org.Hs.eg.db")
  require("clusterProfiler")
  
  enrichment_results <- list()  # Initialize a list to store results for each gene list
  
  for (i in seq_along(gene_lists)) {
    gene_list <- gene_lists[[i]]
    
    # Change ENSEMBL id to ENTREZ
    entrez_gene_list <- mapIds(org.Hs.eg.db, keys = gene_list, keytype = "ENSEMBL", column = "ENTREZID")
    
    # Perform gene enrichment analysis using clusterProfiler
    enrich_result <- enrichKEGG(
      gene          = entrez_gene_list,
      organism      = "hsa",
      pvalueCutoff  = pvalue_cutoff,
      qvalueCutoff  = qvalue_cutoff
    )
    
    # Convert gene ids to gene symbols
    enrich_result2 <- setReadable(enrich_result, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    
    # Store the results in the list
    enrichment_results[[i]] <- enrich_result2
  }
  
  return(enrichment_results)
}





# Cluster Profiler plots together and save them
plot_cluster_profiler <- function(enrich_result, showCategory, plot_title, filename){
  p1 <- as.ggplot(barplot(enrich_result, showCategory=showCategory, title=plot_title))
  p2 <- as.ggplot(dotplot(enrich_result, showCategory=showCategory, title = plot_title))
  p3 <- as.ggplot((cnetplot(enrich_result, showCategory=showCategory, circular=FALSE, title=plot_title)))
  #Concatenate them
  plot_list = list(p1,p2,p3)
  # Arrange and combine the plots using gridExtra
  combined_plots <- grid.arrange(grobs = plot_list, nrow = 2, ncol = 2)
  # Save file
  ggsave(filename, combined_plots,  width = 12, height = 8, dpi = 400)
  print(combined_plots)
}


# Plot UMAP 2D function
create_umap_plot <- function(umap_data, Patient_Colors, plot_title, group, filename) {
  # Create a scatter plot of the UMAP results
  umap_plot <- ggplot(umap_data, aes(x = V1, y = V2, color = umap_data[[group]])) +
    geom_point(size = 6) +
    scale_color_manual(values = Patient_Colors) +
    labs(color = group) +
    xlab("UMAP 1") +
    ylab("UMAP 2") +
    theme_minimal() +
    theme(plot.title = element_text(size = 18), axis.text = element_text(size = 10), legend.text = element_text(size = 10))
  
  # Set the title
  umap_plot <- umap_plot + ggtitle(plot_title)
  
  print(umap_plot)
  
  ggsave(filename, umap_plot, width = 12, height = 8, dpi = 400)
}


# Create a UMAP 3D function
create_UMAP3D <- function(umap_data, Patient_Colors, group) {
  fig <- plot_ly(umap_data, x = ~V1, y = ~V2, z = ~V3, color = umap_data[[group]], colors = Patient_Colors) %>%
    add_markers() %>%
    layout(scene = list(xaxis = list(title = '0'), 
                        yaxis = list(title = '1'), 
                        zaxis = list(title = '2')))
  
  print(fig)
}



# PCA plot function

create_PCA_plot <- function(pca_data, Patient_Colors, plot_title, group, filename, xlab, ylab) {
  # Create a scatter plot of the UMAP results
  pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = pca_data[[group]])) +
    geom_point(size = 6) +
    scale_color_manual(values = Patient_Colors) +
    labs(color = group) +
    xlab(xlab) +
    ylab(ylab) +
    theme_minimal() +
    theme(plot.title = element_text(size = 18), axis.text = element_text(size = 10), legend.text = element_text(size = 10))
  
  # Set the title
  pca_plot <- pca_plot + ggtitle(plot_title)
  
  print(pca_plot)
  
  ggsave(filename, pca_plot, width = 12, height = 8, dpi = 400)
}



# PLot a simple Venn Diagram
#library(VennDiagram)

create_venn_diagram <- function(gene_lists, list_names, output_file) {
  num_lists <- length(gene_lists)
  
  if (num_lists < 2) {
    stop("At least two gene lists are required.")
  }
  
  # Create a color palette for the sets
  #color_palette <- c("red", "green", "blue", "orange", "purple", "pink", "brown", "gray", "yellow", "cyan")
  
  # Generate a Venn diagram with default colors
  venn.plot <- venn.diagram(
    x = gene_lists,
    category.names = list_names,
    filename = output_file,
    output = TRUE,
    #category.colors = color_palette[1:num_lists]
  )
  
  # Customize the colors (fill and border) for each set
  venn.plot <- venn.plot + theme_minimal()
  print(venn.plot)
  
  # Display the Venn diagram
  grid.draw(venn.plot)
}





