
###########################################
## WGCNA on liver transcriptiomics       ##
###########################################
#source("code/vars.R")
#source("code/treatmentEffect.R")


library(WGCNA)
library(DESeq2)
### The following setting is very important
options(stringsAsFactors = FALSE);



##################################################################
#         Read the gene counts table and plot the sample tree    #
##################################################################

###Load file
datExpr<-read.csv("Govaere_2rlog_visualization.csv", header = TRUE, row.names = 1)
# Load metadata
govaere_meta<-read.csv("govaere_meta.csv")




datExpr <- t(datExpr)##Genes as columns
#######################################################
#         Choose soft threshold parameter             #
#######################################################

setwd("C:/MechPath/Masterfile/WGCNA analysis with only common genes/Govaere/October 2022")
Rs = 0.9
powers = c(c(1:20), seq(from=22, to=50, by=2))


sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, RsquaredCut = Rs)
#Scale-free topology fit index as a function of the soft-thresholding power
#pdf(file="sft.pdf", width=9, height=5);
par(mfrow=c(1,2));
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed R^2", ylim=c(-1,1))
     main = paste("Scale independence"));
abline(h=0.9, col="blue")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers, cex=cex1, col="red");

##WGCNA::plotSoftThreshold(sft, Rs, file=paste(myOutPath, "softThresholdGovaere.pdf", sep = ""))




softPower = 10
minModuleSize = 30##30 genes minimum per cluster



#####################################################################
#          Turn data expression into topological overlap matrix     #
#####################################################################

#power = sft$powerEstimate
softPower = 10
adjacency = adjacency(datExpr, power = softPower, type="unsigned")
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
# TOM = TOMsimilarityFromExpr(datExpr, power = power)
dissTOM = 1 -TOM
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method="average");
pdf(file = "gene_cluster.pdf", width=12, height=9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity, log2 Govaere",
     labels=FALSE, hang = 0.04);
dev.off()







#===============================================================================
#
#  Construct modules using adaptive branch pruning of dendrogram
#
#===============================================================================
# Module identification using dynamic tree cut
# We like large modules (clusters), so we set the minimum module size relatively high:
# minModuleSize = 20
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                            pamRespectsDendro = FALSE, minClusterSize = 30)

table(dynamicMods)
length(table(dynamicMods))
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
length(table(dynamicColors))
# Plot the dendrogram and colors underneath
pdf(file="Govaere, 44 clusters, logdata.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05, main= "Gene Dendrogram with 41 modules, Govaere ")

dev.off()




#===============================================================================
#
#  Merge modules
#
#===============================================================================

# Merge close modules based on correlation of eigen genes from modules
MEDissThres=0.20
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
#merge = mergeCloseModules(datExpr, dynamicColors, corFnc = cor, corOptions = list(use='p'), verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
pdf(file = "Govaere merged_Module_Tree_15_clusters.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()



write.table(merge$oldMEs,file="oldMEs.txt");
write.table(merge$newMEs,file="newMEs.txt");




#===============================================================================
#
## Calculate module membership score for each gene (correlation to module eigengene)
#
#===============================================================================
table(merge$colors)
mergedColors = merge$colors
moduleColors = mergedColors

modMembership = as.matrix(cor(datExpr, merge$newMEs, use = 'p'))
modMembershipP = as.matrix(corPvalueStudent(modMembership, nrow(datExpr)))

##Create the starting data frame
geneInfo0 = data.frame(moduleColor = moduleColors,
                       modMembershipP)

write.csv(geneInfo0, "govaere_15 clusters_logdata__30 genes_per cluster.csv")

#===============================================================================
#
#  Plot the heatmap of module eigen-genes (MEs) and samples
#
#===============================================================================

library("pheatmap")

# Heatmap of old module eigen-genes and samples
pdf(file="oldMEs, 44 genes cluster and 0.8 corre.pdf",heigh=80,width=20)
row.names(merge$oldMEs)=names(dataset)
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=6)
dev.off()
# Heatmap of new module eigen-genes and samples
pdf(file="newMEs,  15 genes cluster and 0.8 corre.pdf",heigh=60,width=20)
row.names(merge$newMEs)=names(dataset)
pheatmap(merge$newMEs,cluster_col=T,cluster_row=T,show_rownames=T,show_colnames=T,fontsize=6)
dev.off()



#######################################################################################
#######################################################################################
hubs = chooseTopHubInEachModule(datExpr = datExpr, colorh = moduleColors, omitColors = "grey", power = 9, type = "signed")
tabla = data.frame(hubs)
write.csv(tabla, "hubs_genes_govaere.csv")
#### Get gene symbols from hub genes
annot <-read.csv(file = "hubs_genes_govaere.csv", header = TRUE, row.names = 1)
#######################################################################################
x<-datExpr[,annot$hubs]
x
annot
colnames(x)<-annot$Gene_symbol
top_hubs_govaere<-write.csv(x, "Hub genes for heatmap.csv")


