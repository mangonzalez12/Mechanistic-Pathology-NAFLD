library(umap)
library(ggplot2)
library(gridExtra)
library(cowplot)
library(dplyr)
library(tibble)
library(DESeq2)
library(caret)
library(pheatmap)
library(colorspace)
library(VennDiagram)
library(readxl)
library(RColorBrewer)
library(colorspace)
library(umap)
library(gridExtra)
library(ggplotify)
library(patchwork)
library(circlize)
library(ComplexHeatmap)
list.files()
## Load files 


### Filter whole dataset to extract the hub genes matrix
cts<-read.csv("Govaere_2rlog_visualization.csv", row.names = 1, header = TRUE)
### Get subgroups in a dataframe for later analysis
subgroups <- read.csv("subgroups 6 WGCNA r2log transformation.csv", row.names = 1)
meta <- read.csv("govaere_meta.csv", row.names = 1)
keep<-c("nas_score", "Fibrosis_label")#, "Subgroup")#, "Sample", "Cluster_assignment", "cluster")
meta<- meta[,keep]
## Annotation selected genes (ENS id)
annot <- read.csv("hubs_genes_govaere.csv", row.names = 1, header = TRUE)
x <- as.character(annot$hubs)
df <- cts[x,]
dim(df)
df
# Choose gene names for the plot
select_genes <- annot$Hub_gene
select_genes
df <- t(df)
#colnames(df)==annot$hubs
#colnames(df) <- select_genes
#df
#write.csv(df, "Hub genes for heatmap (clustering)_Govaere from WGCNA output.csv")
#df
# Add subgroups to metafile for heatmap 
rownames(subgroups)==rownames(meta)
meta$Patient_subgroups <- subgroups$x
meta

#################################################################################################
# Genes for heatmap using 14 hub genes
#df <- read.csv("Hub genes for heatmap (clustering)_Govaere from WGCNA output.csv", row.names = 1)
#rownames(meta)==rownames(subgroups)



### 5 levels Fibrosis score, 9 levels NAS-score colors and 6 patient cluster colors
ann_colors <-list(
  nas_score = c(nas_score_0 = "#F4FAFE", nas_score_1 = "#DEEEF7", nas_score_2 = "#C1DBEC", nas_score_3 = "#A1C4E0", nas_score_4 = "#7FABD3", nas_score_5 = "#5C90C6" , nas_score_6 = "#3573B9",  nas_score_7 = "#305596",  nas_score_8 = "#273871"),
  Fibrosis_label = c(Fibrosis_0= "#FADDC3",Fibrosis_1= "#F9C29C", Fibrosis_2= "#F6A173", Fibrosis_3= "#F17B51", Fibrosis_4= "#EA4C3B"),
  Patient_subgroups = c(Subgroup_1= "#023FA5", Subgroup_2="#8C94BF", Subgroup_3= "#D2D3DC", Subgroup_4= "#DDD0D2", Subgroup_5= "#C18692", Subgroup_6="#8E063B")
)
# Colors
nas_score = c(nas_score_0 = "#F4FAFE", nas_score_1 = "#DEEEF7", nas_score_2 = "#C1DBEC", nas_score_3 = "#A1C4E0", nas_score_4 = "#7FABD3", nas_score_5 = "#5C90C6" , nas_score_6 = "#3573B9",  nas_score_7 = "#305596",  nas_score_8 = "#273871")
Fibrosis_label = c(Fibrosis_0= "#FADDC3",Fibrosis_1= "#F9C29C", Fibrosis_2= "#F6A173", Fibrosis_3= "#F17B51", Fibrosis_4= "#EA4C3B")
Patient_subgroups = c(Subgroup_1= "#023FA5", Subgroup_2="#8C94BF", Subgroup_3= "#D2D3DC", Subgroup_4= "#DDD0D2", Subgroup_5= "#C18692", Subgroup_6="#8E063B")


df <- t(df)
dim(df)
## Heatmap
ht <-pheatmap::pheatmap(
  df,
  cluster_rows=T,
  show_rownames=T,
  show_colnames = F,
  cluster_cols=T,
  height = 50, 
  width = 100, 
  border_color = "grey60", 
  clustering_distance_rows = "euclidean", 
  clustering_distance_cols = "euclidean", 
  annotation = meta,
  annotation_colors = ann_colors, ## ann_colors = defined for pathology scores above
  labels_row = select_genes, 
  fontsize_row = 15,
  cutree_cols = 6, 
  fontsize_col = 8, 
  fontsize = 12)

ggsave(ht, filename = "Govaere Heatmap2.png", width = 12, height = 8, dpi = 400)



### Cut heatmap to get the patients in each of the subgroups
#subgroups<- cutree(ht$tree_col, k=6)
#subgroups <- data.frame(subgroups)
#head(subgroups)
#table(subgroups)
### Plot the subgroups
#x11(width = 150, height = 40)
#plot(ht$tree_col, xlab = meta$subgroups)






#######################################################################################
# Plot in a heatmap the pathways/upstream regulators in all the gene clusters (1-14)

canon <- read.csv("02b2.IPAout.CanoP.Gene modules for IPA_Govaere_14 modules.csv", header = TRUE, row.names = 1)
upstream <- read.csv("02b3.IPAout.UpstrP.Gene modules for IPA_Govaere_14 modules_500.csv", header = TRUE, row.names = 1)
dim(canon)

canon <- top_pathways(canon, 50)
upstream <- top_pathways(upstream, 50)
rownames(canon)
# Take top 50 pathways from canon
#canon <- head(canon, 50)
# Take top 50 upstream regulators from upstream
#upstream <- head(upstream,50)
# Scale the canon dataframe with range (0-1)
preproc <- preProcess(canon, method = c("range"))
sc_canon <- predict(preproc, newdata = canon)

# Scale the upstream dataframe with range (0-1)
preproc1 <- preProcess(upstream, method = c("range"))
sc_upstream <- predict(preproc1, newdata = upstream)

# Color
pal <- colorRampPalette(c("white", "blue"))
palette <- pal(10)


# Plot size canon file 12 x 13, respectively is good
width_canon = 12
height_canon =13

#Top100
#width_canon = 15
#height_canon =40

#Top200
#width_canon = 15
#height_canon =80

#colnames_clusters <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "Cluster 11", "Cluster 12", "Cluster 13", "Cluster 14")
colnames_clusters <- c("Module 1", "Module 2", "Module 3", "Module 4", "Module 5", "Module 6", "Module 7", "Module 8", "Module 9", "Module 10", "Module 11", "Module 12", "Module 13", "Module 14")
# Reset colnames of both dataframes with gene module names
colnames(sc_canon) <- colnames_clusters
colnames(sc_upstream) <- colnames_clusters


# Heatmap Canon
ht_canon <-pheatmap::pheatmap(
  sc_canon,
  cluster_rows=F,
  show_rownames=T,
  cluster_cols=F,
  border_color = "grey60", 
  color = colorRampPalette(c("white", "navy blue"))(20), ## ann_colors = defined for pathology scores above
  #fontsize_row = 16,
  #fontsize_col = 16, 
  #fontsize = 16,
  #cellwidth = 16,
  #cellheight = 16,  
  fontsize_row = 11,
  fontsize_col = 11, 
  fontsize = 12,
  cellwidth = 10,
  cellheight = 10,
  width = width_canon,
  height = height_canon,
  legend = FALSE,
  main = "A)     Canonical Pathways in gene modules"
  )

pdf("Canonical_Pathways_IPA_gene_modules_Feb22.pdf", width = width_canon, height = height_canon)
print(ht_canon)
dev.off()

print(ht_canon)

# Save a PNG plot
ggsave("Canonical_Pathways_IPA_gene_modules_Feb22.png", plot = ht_canon, width = width_canon, height = height_canon, dpi=400)


# Plot size upstream file
width_upstream = 6
height_upstream =13

# Convert dataframe to matrix
sc_upstream <- as.matrix(sc_upstream)
sc_upstream
# Heatmap upstream
ht_upstream <-pheatmap::pheatmap(
  sc_upstream,
  cluster_rows=F,
  show_rownames=T,
  cluster_cols=F,
  border_color = "grey60", 
  color = colorRampPalette(c("white", "navy blue"))(20), ## ann_colors = defined for pathology scores above
  fontsize_row = 11,
  fontsize_col = 11, 
  fontsize = 12,
  cellwidth = 10,
  cellheight = 10,
  width = width_upstream,
  height = height_upstream,
  legend = TRUE,
  main = "B)     Upstream regulators in gene modules"
)
# Plot upstream
print(ht_upstream)

# Save upstream as PDF
pdf("Upstream_IPA_gene_modules_Feb22.pdf", width = width_upstream, height = height_upstream)
print(ht_upstream)
dev.off()

# Save a PNG plot
ggsave("Upstream_IPA_gene_modules_Feb22.png", plot = ht_upstream, width = width_upstream, height = height_upstream, dpi=400)

# Arrange the two heatmaps 
# Save plots as ggplots with ggplotify package
p1 <- as.ggplot(ht_canon, hjust = 0, vjust = 0)
p2 <- as.ggplot(ht_upstream, hjust = 0, vjust = 0)

# Use patchwork to put them together
combined_plot <- p1 + theme(plot.margin = unit(c(0,192,0,0), "pt")) + p2 + 
  plot_layout(widths = c(4, -1.9,4.5),guides = "collect")& theme(legend.position = "top")
# Print plot
print(combined_plot)


width= 20
height=12

# Save upstream as PDF
pdf("Combined plot canon and upstream 14 gene modules_Feb22.pdf", width = width, height = height)
print(combined_plot)
dev.off()
# Save a PNG plot
ggsave("Combined plot canon and upstream 14 gene modules_Feb22.png", plot = combined_plot, width = width, height = height, dpi=1000)







#########################################################################################
#########################################################################################
## Plot in a heatmap the pathways/upstream regulators in all the 6 Patient subgroups
canon <- read.csv("04d2.Cano.Zclustered.6comps_manuallycurated.csv", header = TRUE, row.names = 1)
upstream <- read.csv("04d4.UpstrDrugs.Zclustered.6comps.csv", header = TRUE)
upstream <- read.csv("04d4.UpstrDrugs.Zclustered.6comps_Top50.csv", header = TRUE, row.names = 1)

#Filter low count pathways
canon <- remove_rows_absolute_below_threshold(canon, threshold = 2.5, min_columns =4)
canon
dim(canon)
dim(upstream)
# Take top 50 pathways from canon
canon <- head(canon, 50)
canon <- head(canon, 100)
canon
# Color for blue to red
pal<-choose_palette()
Diverging_Colors <- pal(9)
# Color
Diverging_Colors <- rev(brewer.pal(9, "RdBu"))

# Colors used
Diverging_Colors <- c("#00669C", "#0090B3", "#AAC6D2", "#DBB9C1", "#C26680", "#9F0045")

# Plot size canon file 12 x 13, respectively is good
#All pathways
width_canon = 12
height_canon =80
layout_matrix <- matrix(c(2, 1), ncol = 2, byrow = TRUE)
layout_heights <- c(4, 1)
#Top50
width_canon = 12
height_canon =13
#Top100
width_canon = 15
height_canon =40

#Check column names
colnames(canon)
# Colnames change
colnames_patient_subgroups <- c("Subgroup 1", "Subgroup 2", "Subgroup 3", "Subgroup 4", "Subgroup 5", "Subgroup 6")
# Reset colnames of both dataframes with patient subgroup names
colnames(canon) <- colnames_patient_subgroups
colnames(upstream) <- colnames_patient_subgroups
canon
dim(canon)
#upstream





## Cut the plot in 2 for the figure in the manuscript
# Set up the layout with two rows
layout(matrix(1:2, nrow = 2))


# Heatmap Canon 1st part
ht_canon1 <-pheatmap(
  canon[1:floor(nrow(canon)/2),],
  cluster_rows=F,
  show_rownames=T,
  cluster_cols=F,
  border_color = "grey60", 
  color = Diverging_Colors, ## ann_colors = defined for pathology scores above
  fontsize_row = 16,
  fontsize_col = 16,
  fontsize = 16,
  cellwidth = 25,
  cellheight = 25,
  #fontsize_row = 9,
  #fontsize_col = 9, 
  #fontsize = 12,
  #cellwidth = 8,
  #cellheight = 8,
  width = width_canon,
  height = height_canon,
  legend = FALSE,
  display_numbers = TRUE,
  fontsize_number = 9,
  main = "    A)  Canonical Pathways in Patient Subgroups part 1",
  mar = c(10, 5, 5, 15)  # Adjust the bottom margin (mar[4]) to make more space for row labels
)


pdf("Canonical Pathways in Patient Subgroups part 1.pdf", width = width_canon, height = height_canon)
print(ht_canon1)
dev.off()


print(ht_canon1)

# Save a PNG plot
ggsave("Canonical Pathways in Patient Subgroups part 1.png", plot = ht_canon1, width = width_canon, height = height_canon, dpi=400, limitsize = FALSE)




# Heatmap Canon 2nd part
ht_canon2 <-pheatmap(
  canon[((nrow(canon)/2) + 1):nrow(canon), ],
  cluster_rows=F,
  show_rownames=T,
  cluster_cols=F,
  border_color = "grey60", 
  color = Diverging_Colors, ## ann_colors = defined for pathology scores above
  fontsize_row = 16,
  fontsize_col = 16,
  fontsize = 16,
  cellwidth = 25,
  cellheight = 25,
  #fontsize_row = 9,
  #fontsize_col = 9, 
  #fontsize = 12,
  #cellwidth = 8,
  #cellheight = 8,
  width = width_canon,
  height = height_canon,
  legend = TRUE,
  display_numbers = TRUE,
  fontsize_number = 9,
  main = "              B)  Canonical Pathways in Patient Subgroups part 2",
  mar = c(10, 5, 5, 15)  # Adjust the bottom margin (mar[4]) to make more space for row labels
)

# Combine Canonical part 1 and part 2
p1 <- as.ggplot(ht_canon1, hjust = 0, vjust = 0)
p2 <- as.ggplot(ht_canon2, hjust = 0, vjust = 0)

# Use patchwork to put them together
combined_plot <- p1 +
  plot_spacer() + 
  p2 + 
  plot_layout(widths = c(4, -1.9,4.5),guides = "collect")+
  theme(legend.position = "top")
# Print plot
print(combined_plot)

# Set size
width= 23
height=20

# Save upstream as PDF
pdf("Combined plots in Patient Subgroups_Figure paper_manually curated.pdf", width = width, height = height)
print(combined_plot)
dev.off()


# Save a PNG plot
ggsave("Combined plots in Patient Subgroups_Figure paper_manually curated.png", plot = combined_plot, width = width, height = height, dpi=1000)







# Size PsP  
#Top50
width_upstream = 12
height_upstream =13

upstream
# Heatmap upstream
ht_upstream <-pheatmap(
  upstream,
  cluster_rows=F,
  show_rownames=T,
  cluster_cols=F,
  border_color = "grey60", 
  color = Diverging_Colors, ## ann_colors = defined for pathology scores above
  fontsize_row = 16,
  fontsize_col = 16, 
  fontsize = 16,
  cellwidth = 20,
  cellheight = 20,
  width = width_upstream,
  height = height_upstream,
  legend = TRUE,
  main = "Upstream regulators in Patient Subgroups Top50"
)
# Plot upstream
print(upstream)
dim(upstream)
# Save upstream as PDF
pdf("Upstream IPA regulators in Patient Subgroups_Top50.pdf", width = width_upstream, height = height_upstream)
print(ht_upstream)
dev.off()

# Save a PNG plot
ggsave("Upstream IPA regulators in Patient Subgroups_Top50.png", plot = ht_upstream, width = width_upstream, height = height_upstream, dpi=400)


# Arrange the two heatmaps 
# Save plots as ggplots with ggplotify package
p1 <- as.ggplot(ht_canon, hjust = 0, vjust = 0)
p2 <- as.ggplot(ht_upstream, hjust = 0, vjust = 0)

# Use patchwork to put them together
combined_plot <- p1 + plot_spacer() + p2 + 
  plot_layout(widths = c(4, -1.9,4.5),guides = "collect")& theme(legend.position = "top")
# Print plot
print(combined_plot)


width= 25
height=20

# Save upstream as PDF
pdf("Combined plots in Patient Subgroups.pdf", width = width, height = height)
print(combined_plot)
dev.off()


# Save a PNG plot
ggsave("Combined plots in Patient Subgroups.png", plot = combined_plot, width = width, height = height, dpi=1000)










#############################################################################
#############################################################################
#############################################################################
#############################################################################

#### CIRCULAR HEATMAPS



## Read files canon and zscore from DESEQ analysis one vs Rest
canon <- read.csv("04d1.Cano.Pclustered.6comps.csv", row.names = 1, header = TRUE, sep = ",")
zscore <- read.csv("04d2.Cano.Zclustered.6comps.csv", row.names = 1, header = TRUE, sep = ",")
# Filter both canon and zscore by common pathways in both
common_pathways <- intersect(rownames(canon), rownames(zscore))
canon <- canon[common_pathways, , drop = FALSE]
zscore <- zscore[common_pathways, , drop = FALSE]
# Read File with pathways categories
szscore <- read.csv("IPA_zscore_all_split_A_D.csv", row.names = 1, header = TRUE, sep = ",")
# Filter both canon and zscore with categories from A to D
relevant_pathways <- intersect(rownames(szscore), rownames(zscore))
canon1 <- canon[relevant_pathways, , drop = FALSE]
zscore1 <- zscore[relevant_pathways, , drop = FALSE]
#Check dimensions
dim(canon)
dim(zscore)
dim(canon1)
dim(zscore1)
rownames(canon1)==rownames(zscore1)
# Filter pathways
selected_pathways <- filter_pathways(canon1, zscore1, threshold_significance = 90, threshold_directionality = 0.5)
length(selected_pathways)

canon2 <- canon[selected_pathways,]
zscore2 <- zscore[selected_pathways,]
szscore2 <- szscore[selected_pathways,]
table(szscore2$Group)
zscore2

# Ensure splitz2$Group is a factor
szscore2$Group <- factor(szscore2$Group)
# Assign custom colors to factor levels
vec_col <- szscore2$Group
vec_col
table(vec_col)
col_fun = colorRamp2(c(-4, 0, 4), c("blue", "white", "red"))
col_mechanisms <- list("A" = "A: Fibrosis",
                       "B" = "B: Lipid Metabolism",
                       "C" = "C: Inflammation",
                       "D" = "D: Glucose metabolism")

# Define colors for each sector
sector_colors <- c("black", "deepskyblue3", "darkblue", "darkorange3")
# Repeat each color for the number of rows in each sector
names(sector_colors) <- unique(vec_col)  # Assign colors to unique classes
# Assign colors to row names based on vec_col
sector_row_colors <- sector_colors[vec_col]
# Change column names
colnames <- c("Sub 1", "Sub 2", "Sub 3", "Sub 4", "Sub 5", "Sub 6")
colnames(zscore2) <- colnames
## Clear
circos.clear()
graphics.off()

pdf("Ingenuity Pathway Analysis on 6 patient Subgroups.pdf", height = 9, width = 12)

# Define parameters
start_degree <- 90
gap_degree <- c(4, 4, 4, 8)
track_margin <- c(0, 0)
legend_title_directionality <- "Directionality"
legend_title_mechanisms <- "Mechanisms"

# Circular Heatmap
circos.par(
  start.degree = start_degree,
  gap.degree = gap_degree,
  track.margin = track_margin
)

# Legend
legend_directionality <- Legend(
  title = legend_title_directionality,
  col_fun = col_fun,
  grid_height = unit(4, "mm"),
  grid_width = unit(4, "mm")
)
legend_mechanisms <- Legend(
  title = legend_title_mechanisms,
  col_mechanisms,
  grid_height = unit(8, "mm"),
  grid_width = unit(8, "mm")
)

# Set row names
rownames(zscore2) <- factor(rownames(zscore2))

# Circular Heatmap
circos_ht <- circos.heatmap(
  zscore2,
  split = vec_col,
  col = col_fun,
  cluster = TRUE,
  rownames.side = "inside",
  cell.border = "black",
  bg.lwd = 1,
  bg.lty = 1,
  cell.lwd = 1,
  track.height = 0.15,
  show.sector.labels = TRUE,
  dend.side = "outside",
  dend.track.height = 0.1,
  rownames.cex = 0.5,
  cell.lty = 1,
  bg.border = "gray50",
  rownames.col = sector_row_colors,
  rownames.font = par("font")
)

# Adding a new track for column labels
circos.track(
  track.index = get.current.track.index() + 1 ,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 4) { # the last sector
      column_names <- colnames(zscore2)
      num_columns <- length(column_names)
      ## Set the order correctly of the columns (Patient Subgroups)
      column_names <- rev(column_names)
      
      circos.text(
        rep(CELL_META$cell.xlim[2], num_columns) + convert_x(1, "mm"),
        1:num_columns - 0.5,
        column_names,
        cex = 0.7,
        adj = c(0, 0.8),
        facing = "inside"
      )
    }
  },
  bg.border = NA
)

# Draw legends
draw(legend_directionality, x = unit(0.80, "npc"), y = unit(0.95, "npc"), just = c("right", "top"))
draw(legend_mechanisms, x = unit(0.20, "npc"), y = unit(0.80, "npc"), just = c("right", "top"))

# Add title
title("Ingenuity Pathway Analysis on DEGs from 6 Patient Subgroups")

dev.off()








## 14 gene clusters circular heatmap
# Plot in a heatmap the pathways/upstream regulators in all the gene clusters (1-14)

canon <- read.csv("02b2.IPAout.CanoP.Gene modules for IPA_Govaere_14 modules.csv", header = TRUE, row.names = 1)
upstream <- read.csv("02b3.IPAout.UpstrP.Gene modules for IPA_Govaere_14 modules_500.csv", header = TRUE, row.names = 1)

# Rename columns
colnames_clusters <- c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9", "Cluster 10", "Cluster 11", "Cluster 12", "Cluster 13", "Cluster 14")
# Reset colnames of both dataframes with gene module names
colnames(canon) <- colnames_clusters
colnames(upstream) <- colnames_clusters

# Filter canon based on szscore file with pathway categories
x <- intersect(rownames(canon), rownames(szscore))
canon <- canon[x,]

# See significance of those pathways and select top XX number of pathways
canon2 <- top_pathways(canon, threshold_significance = 100)
canon2

# Filter szscore to have the same as the top pathways
szscore2 <- szscore[rownames(canon2),]


# Ensure splitz2$Group is a factor
szscore2$Group <- factor(szscore2$Group)
# Assign custom colors to factor levels
vec_col <- szscore2$Group
vec_col
table(vec_col)
col_fun = colorRamp2(c(10, 0), c("blue", "white"))
legend_title_significance <- "Significance"
col_mechanisms <- list("A" = "A: Fibrosis",
                       "B" = "B: Lipid Metabolism",
                       "C" = "C: Inflammation",
                       "D" = "D: Glucose metabolism")
significance <- range(canon2)
significance
## Clear
circos.clear()
graphics.off()

pdf("Ingenuity Pathway Analysis on 14 gene modules100v3.pdf", height = 12, width = 15, pointsize = 12)

# Define parameters
start_degree <- 90
gap_degree <- c(4, 4, 4, 12)
track_margin <- c(0, 0)
legend_title_mechanisms <- "Mechanisms"

# Circular Heatmap
circos.par(
  start.degree = start_degree,
  gap.degree = gap_degree,
  track.margin = track_margin
)

# Legend
legend_significance <- Legend(
  title = legend_title_significance,
  col_fun = col_fun,
  grid_height = unit(4, "mm"),
  grid_width = unit(4, "mm")
)
legend_mechanisms <- Legend(
  title = legend_title_mechanisms,
  col_mechanisms,
  grid_height = unit(8, "mm"),
  grid_width = unit(8, "mm")
)

# Set row names
rownames(canon2) <- factor(rownames(canon2))

# Circular Heatmap
circos_ht <- circos.heatmap(
  canon2,
  split = vec_col,
  col = col_fun,
  cluster = TRUE,
  rownames.side = "inside",
  cell.border = "black",
  bg.lwd = 1,
  bg.lty = 1,
  cell.lwd = 1,
  track.height = 0.15,
  show.sector.labels = TRUE,
  dend.side = "outside",
  dend.track.height = 0.1,
  rownames.cex = 0.7,
  cell.lty = 1,
  bg.border = "gray50",
  rownames.col = 1:nrow(canon2) %% 10 + 5,
  rownames.font = par("font")
)

# Adding a new track for column labels
circos.track(
  track.index = get.current.track.index() + 1 ,
  panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 4) { # the last sector
      column_names <- colnames(canon2)
      num_columns <- length(column_names)
      ## Set the order correctly of the columns (Patient Subgroups)
      column_names <- rev(column_names)
      
      circos.text(
        rep(CELL_META$cell.xlim[2], num_columns) + convert_x(1, "mm"),
        1:num_columns - 0.5,
        column_names,
        cex = 0.4, #0.7 top100, #0.4 top50
        adj = c(0, 1),
        facing = "inside"
      )
    }
  },
  bg.border = NA
)

# Draw legends
draw(legend_mechanisms, x = unit(0.18, "npc"), y = unit(0.80, "npc"), just = c("right", "top"))
draw(legend_significance, x = unit(0.85, "npc"), y = unit(0.85, "npc"), just = "right")

# Add title
title("Ingenuity Pathway Analysis on 14 gene modules")

dev.off()





