### Gene Set Enrichment analysis using cluster profiler
library(clusterProfiler)
library(fgsea)
library(org.Hs.eg.db)
library(ComplexUpset)
library(readxl)
library(fgsea)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(utils)
library(readr)
library(dplyr)
library(pheatmap)
library(biomaRt)
library(igraph)
library(openxlsx)


#####Functions for all analyses

## A function to round cts
round_df <- function(x, digits){
  numeric_columns<- sapply(x, mode)=='numeric'
  x[numeric_columns] <- round(x[numeric_columns], digits)
}

## A function to remove decimals
remove_decimals <- function(cts){
  genes_with_decimals <-rownames(cts)
  genes_without_decimals<-gsub("\\..*", "", genes_with_decimals)
  rownames(cts) <- genes_without_decimals
  return(cts)
}


### Make gene symbol database from ENSEMBLS in study
makeGeneSymbols_db <- function(data_frame) {
  require("biomaRt")
  
  # Get gene Symbols in all study
  ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  attributes <- c("ensembl_gene_id", "external_gene_name")

  # Retrieve gene symbols based on Ensembl IDs from your data frame
  gene_ids <- rownames(data_frame)
  gene_info <- getBM(attributes, filters = "ensembl_gene_id", values = gene_ids, mart = ensembl)
  colnames(gene_info) <- c("ENSEMBL_ID", "Gene_Symbol")
  
  return(gene_info)
}


#Get gene SYMBOLS
map_gene_symbols_from_db <- function(ensembl_ids, gene_data) {
  match_indices <- match(rownames(ensembl_ids), gene_data$ENSEMBL_ID)
  gene_symbols <- gene_data$Gene_Symbol[match_indices]
  return(gene_symbols)
}

#Convert and replace
convertAndReplace <- function(list_df) {
  for (i in seq_along(list_df)) {
    char_cols <- sapply(list_df[[i]], is.character)
    char_cols[length(char_cols)] <- FALSE
    
    list_df[[i]][, char_cols] <- lapply(list_df[[i]][, char_cols, drop = FALSE], function(x) as.numeric(gsub(",", ".", x)))
  }
  return(list_df)
}


# Function to extract Symbols and Logfoldchange
get_symbol_logFC <- function(df) {
  columns <- c("Symbol", "log2FoldChange")
  df <- df[,columns]
  gene_list = df$log2FoldChange
  names(gene_list) = df$Symbol
  gene_list = sort(gene_list, decreasing = TRUE)
  gene_list = gene_list[!duplicated(names(gene_list))]
  gene_list2 <- gene_list[names(gene_list) != ""]
  return(gene_list2)
}


# Use Apply the function to each dataframe in the 
#list = new_dataframes <- lapply(list_df, get_symbol_logFC )

##### Functions for data ordering
##################################################################################################################################
##################################################################################################################################
### Write a function to cut a tree and obtain a dataframe with the sample subgroups
get_subgroups<- function(tree_hc, k){
  ####tree_hc<-hclust() object, either on columns or rows
  groups <- cutree(tree_hc, k)
  return(groups)
}

### Write a function to subset the gene expression matrix
subset_cts <- function(cts, vector){
  cts.new <- cts[vector,]
  return(cts.new)
}

### Write a function to order dataframe based on a column
order_rows<-function(meta, order_vector){### input data metafile and vector to order the rownames
  df_ordered<- meta[order(order_vector),]### Sort metadata rownames (samples) based on vector F score
  return(df_ordered)### Return meta ordered
}

### Write a function to order cts file using the meta_ordered file rownames
order_cts <- function(cts, order_vector){
  ordered_cts <- cts[order_vector]
  return(ordered_cts)
}

### Function to round values 
round_df <- function(x, digits) {
  # round all numeric variables
  # x: data frame 
  # digits: number of digits to round
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

####Function to remove decimals from the gene names using regex
remove_decimals <-function(cts){
  genes_with_decimals<-rownames(cts)
  genes_without_decimals<- gsub("\\..*", "", genes_with_decimals)
  rownames(cts)<-genes_without_decimals
  return(cts)
}

###Preprocessing of the data
##################################################################################################################################
##################################################################################################################################
## Function to filter cts file for low-gene count values (mainly zero genes)

###Filter dataset for low-counts
###>75% of the samples has values >=10
###All the gene-rows that have <=75% with values >=10 will be removed.
####I would not remove values <10 but leave them in, as they are.

###I agree with Martien that you should remove all rows (genes) of which
####a percentage (e.g. 75%) of the samples do not have more than 10 counts.

filter_low_counts<-function(cts, filter){
  #filter = 5 or 10
  threshold <- round(length(cts) * 0.75)
  filtered_cts <- cts[rowSums(cts>=filter) >= threshold,]
  return(filtered_cts)
}


### Write a function to normalize data as follows: 
##   nCts = (raw_gene_value/total sum each sample)*average(total sum of samples)
normalize <- function(cts){
  vector <- colSums(cts)
  average_sum<-mean(colSums(cts))
  cts.sweep <- sweep(cts, MARGIN = 2, STATS =vector, '/') 
  cts.norm   <- cts.sweep * average_sum
  return(cts.norm)
}

### Write a function to normalize for healthy samples as follows:
##    2logR = log2(nCts/mean(healthy sample average per row (genes from healthy samples)))

transformation_Rlog_<-function(cts.norm){
  boolean_healthy<-meta$Fibrosis_label=="Fibrosis_0"
  healthy_samples= cts.norm[,boolean_healthy]
  healthy_samples_means <-rowMeans(healthy_samples)
  cts.r2log<-log2((cts.norm+2)/ (healthy_samples_means+2)) ### Add integer to avoid negative values
  return(cts.r2log)
}

transformation_Rlog_positive_values<-function(cts.norm){
  boolean_healthy<-meta$Fibrosis_label=="Fibrosis_0"
  healthy_samples= cts.norm[,boolean_healthy]
  healthy_samples_means <-rowMeans(healthy_samples)
  cts.r2log <- log2((cts.norm+2)/ (healthy_samples_means+2))+3### Add integer to avoid negative values
  return(cts.r2log)
}

#### Write a function to scale the data for dimensionality reduction
## scaling function as following formula:
##  scaled_value = (value - mean value)/standard deviation
scale_cts<-function(cts.r2log){
  df<-t(cts.r2log)
  mean<-mean(colSums(df))##per row genes
  sd<- sd(colSums(df))
  cts.scaled <- cts.r2log - mean/ sd
  return(t(cts.scaled))
}



################# Functions to filter pathway analysis

# Filter IPA output from 14 gene modules

remove_rows_with_zeros <- function(dataframe, threshold = 7) {
  # Count the number of zeros in each row
  zero_counts <- apply(dataframe, 1, function(row) sum(row == 0))
  
  # Identify rows with 7 or more zeros
  rows_to_remove <- which(zero_counts >= threshold)
  
  # Remove identified rows
  new_dataframe <- dataframe[-rows_to_remove, ]
  
  return(new_dataframe)
}


# Function to filter the low values in the heatmap to plot the most relevant
# canonical pathways from IPA analysis

remove_rows_absolute_below_threshold <- function(data, threshold, min_columns) {
  # Function to remove rows where absolute value is below threshold in four or more columns
  
  # Apply the condition to each row
  row_condition <- apply(abs(data) < threshold, 1, function(row) sum(row) >= min_columns)
  
  # Filter rows based on the condition
  filtered_data <- data[!row_condition, ]
  
  return(filtered_data)
}



# GeNeck extract top hubs from each network file
# Create a graph from the data
extract_hubs <- function(network, num_hubs){
  graph <- graph_from_data_frame(network, directed = FALSE)
  # Calculate node degrees
  degrees <- degree(graph)
  # Identify hub genes (nodes with high degree)
  hub_genes <- names(sort(degrees, decreasing = TRUE)[1:num_hubs])
  df <- data.frame(hub_genes)
  # Return df with hubs
  return(df)
}

## Write a excel file with hub gene names
createExcelWorkbook <- function(hub_lists, file_name = "output.xlsx") {
  wb <- createWorkbook()
  
  for (i in seq_along(hub_lists)) {
    sheet_name <- paste0("Net", i)
    addWorksheet(wb, sheetName = sheet_name)
    writeData(wb, sheet = sheet_name, x = hub_lists[[i]], startCol = 1, startRow = 1, colNames = TRUE)
  }
  
  saveWorkbook(wb, file = file_name)
}

filter_pathways <- function(canon_df, zscore_df, threshold_significance = 100, threshold_directionality = 0.5) {
  # Calculate row sums of canon_df to rank pathways
  pathway_ranks <- rowSums(canon_df)
  
  # Select top pathways based on threshold_significance
  top_pathways <- names(head(sort(pathway_ranks, decreasing = TRUE), threshold_significance))
  
  # Find pathways with directionality higher than absolute value of threshold_directionality in zscore_df
  directionality_condition <- apply(abs(zscore_df[top_pathways,]) > threshold_directionality, 1, any, na.rm = TRUE)
  
  # Combine the two conditions to get the final set of pathways
  selected_pathways <- top_pathways[directionality_condition]
  
  return(selected_pathways)
}



# Select top pathways
top_pathways <- function(canon_df, threshold_significance = 100) {
  # Calculate row sums of canon_df to rank pathways
  pathway_ranks <- rowSums(canon_df)
  
  # Select top pathways based on threshold_significance
  top_pathways <- names(head(sort(pathway_ranks, decreasing = TRUE), threshold_significance))
  # Filter canon with top pathways
  canon_filtered <- canon_df[top_pathways,]
  
  return(canon_filtered)
}
