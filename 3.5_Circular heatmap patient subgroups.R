library(circlize)
library(ComplexHeatmap)

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
# Larger font size pathways
selected_pathways <- c("Collagen biosynthesis and modifying enzymes",
                       "IL-6 Signaling",
                       "IL-17 Signaling",
                       "NAFLD Signaling Pathway",
                       "Platelet homeostasis",
                       "Insulin Secretion Signaling Pathway",
                       "AMPK Signaling")

# Define row names cex vector
rownames_cex <- ifelse(rownames(zscore2) %in% selected_pathways, 0.8, 0.53)


## Clear
circos.clear()
graphics.off()

pdf("Ingenuity Pathway Analysis on 6 patient Subgroupsv7.pdf", height = 9, width = 12)

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
  rownames.cex = rownames_cex,
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