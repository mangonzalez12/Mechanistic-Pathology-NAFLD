# preprocessing RNA seq data + ordinal regression

# install packages
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#BiocManager::install("edgeR")
#install.packages('ordinal')
# install.packages('ggeffects')

# load packages
library(edgeR)
library(ordinal)
library(ggeffects)

# load data
govaere<-read.csv("Govaere_2rlog_visualization_6_subgroups for OR.csv", header = TRUE, row.names = 1)
# Load metadata
govaere_meta<-read.csv("govaere_meta.csv", header = TRUE, row.names = 1)

##Non negative values
govaere <- govaere + 10
range(govaere)
# create DGElist
group <- govaere_meta$Cluster_WGCNA
y <- DGEList(counts= t(govaere), group = group)
y$samples
y$counts

# filter out genes with low or no expression
# minimum library size in milions = L = 7 --> 10/7 = 1.4, smallest group has 6 samples (Chen et al., 2016)
keep <- rowSums(cpm(y) > 1.4) >= 6
table(keep) # 17446 genes
y <- y[keep, , keep.lib.sizes=FALSE]
nrow(y$counts)

### normalize counts TMM normalization
y <- calcNormFactors(y, method = 'TMM')
y$samples

# log2(CPM) values = gene-level abundance estimate
logCPM <- cpm(y, normalized.lib.sizes=TRUE, log=TRUE) # filtered normalized gene abundance estimates

#####################################################################################
##################### MDS plot for fibrosis stage ##############################################

# color scheme for fibrosis stage
govaere_meta$Cluster_WGCNA <- as.factor(govaere_meta$Cluster_WGCNA)
levels(govaere_meta$Cluster_WGCNA)
col.cell <- c("purple","orange", "red", "black", "blue", "yellow")[govaere_meta$Cluster_WGCNA]
data.frame(govaere_meta$Cluster_WGCNA,col.cell)
plotMDS(y,labels = govaere_meta$Cluster_WGCNA, pch = 19, col=col.cell)
# Add a title
title("MDS plot with Subgroup in colors")

##########################################################################
#################### ordinal regresssion #################################

str(govaere_meta)
#Covariates$sex <- as.numeric(Covariates$sex)

# describe the data that went into the analysis
#summary(Covariates)      
#table(Covariates$sex)    
#mean(Covariates$age); sd(Covariates$age)
#mean(Covariates$BMI); sd(Covariates$BMI) 
#mean(Covariates$`WC (cm)` , na.rm = T); sd(Covariates$`WC (cm)`, na.rm = T) 
#table(Covariates$T2DM)
#table(Covariates$HT)
#table(Covariates$DL)
#table(Covariates$overall_score)

# check log2(CPM) values = gene-level abundance estimate
str(logCPM)
dim(logCPM)
#rownames(govaere_meta) <- Covariates$FastqID 
table(rownames(govaere_meta)==colnames(logCPM))

govaere_meta$Cluster_WGCNA <- as.ordered(govaere_meta$Cluster_WGCNA)
govaere_meta$GENEi <- rep(NA,times=nrow(govaere_meta))
head(govaere_meta)

# ordinal model: clm(response ~ predictors)

ordinal_model <- function(data)
{
clm(govaere_meta$Cluster_WGCNA ~ GENEi, data = data, link = "logit")
}

# Testrun for 1 gene to determine the number of estimates in output
govaere_meta$GENEi <- logCPM[10,]
ord <- ordinal_model(govaere_meta)
coeff <- summary(ord)$coefficients  
print(coeff)

#######################################################################
######## fit ordinal regression model for all genes ###################
# Create objects to save the output from ordinal regression
Nestimates <-nrow(coeff)
Ngene      <-nrow(logCPM)
pval       <-matrix(NA,Ngene,Nestimates)
estimate   <-matrix(NA,Ngene,Nestimates)
SE         <-matrix(NA,Ngene,Nestimates)
z_value    <-matrix(NA,Ngene,Nestimates)
N          <- matrix(NA,Ngene,1)
Ngene
Nestimates

# Run model for all GENEs and store output
for (i in 1:Ngene)
{
  govaere_meta$GENEi <- logCPM[i,]
  N[i,] <- length(which(!is.na(govaere_meta$GENEi)))
  coeff <- summary(ordinal_model(govaere_meta))$coefficients
  for (j in 1: Nestimates)
  {
    estimate[i,j] <- coeff[j,1]
    SE[i,j] <- coeff[j,2]
    z_value[i,j]  <- coeff[j,3]
    pval[i,j]     <- coeff[j,4]
  }
}

# assign gene names as rownames
rownames(estimate) <-  rownames(logCPM)
rownames(SE) <-  rownames(logCPM)
rownames(z_value) <-  rownames(logCPM)
rownames(pval) <-  rownames(logCPM)
# assign column name
colnames(estimate) <- rownames(coeff)
colnames(SE) <- rownames(coeff)
colnames(z_value) <- rownames(coeff)
colnames(pval) <- rownames(coeff)

# save output ordinal regression
save(estimate,pval,SE,z_value,N, file= 'output ordinal regression') 

##############################################################################
#################### significant genes #######################################

# p-value correction for multiple testing
# bonferoni correction
bonf_p <- 0.05/(Ngene)
# number of significant genes 
table(pval[,6] < bonf_p) ###################################### 850 significant genes

p_sig_genes <- subset(pval, pval[,6] < bonf_p)
estimates_sig_genes <- estimate[row.names(p_sig_genes),]
se_sig_genes <- SE[row.names(p_sig_genes),]
z_value_sig_genes <- z_value[row.names(p_sig_genes),]

sig_genes <- cbind(p_sig_genes[,6], estimates_sig_genes[,6], se_sig_genes[,6], z_value_sig_genes[,6])
colnames(sig_genes) <- c('p-value', 'estimates', 'SE', 'z-value')
save(sig_genes, file = 'significant genes')

genes_normalized <- y$counts
sig_genes_normalized <- genes_normalized[rownames(sig_genes),]
save(sig_genes_normalized, file = 'significant genes normalized counts')


# False discovery threshold (FDR)
q = 0.05
q_values <- p.adjust(pval[,6], method = "hochberg")
table(q_values < q) ######################################### 58 significant genes
pval = cbind(pval, q_values)

q_sig_genes <- subset(pval, pval[,6] < q)
q_val_sig_genes <- q_sig_genes[,-1:-5]
qestimates_sig_genes <- estimate[row.names(q_sig_genes),]
qse_sig_genes <- SE[row.names(q_sig_genes),]
qz_value_sig_genes <- z_value[row.names(q_sig_genes),]

q_sig_genes <- cbind(q_val_sig_genes, qestimates_sig_genes[,6], qse_sig_genes[,6], qz_value_sig_genes[,6])
colnames(q_sig_genes) <- c('p-value', 'sex', 'q_value', 'estimates', 'SE')

########################################################################################

# top 400 genes
load('output ordinal regression.rds')
dim(pval)
colnames(pval)

genes <- cbind(pval[,5], estimate[,5], SE[,5], z_value[,5])
colnames(genes) <- c('pvalue', 'estimates', 'SE', 'z-value')
genes <- as.data.frame(genes)
write.csv(genes, "All genes Govaere OR.csv")

head(genes)
sorted_pval <- genes[order(genes$pvalue),]
top_400 <- sorted_pval[1:400,]
dim(top_400)
write.csv(top_400, "top_400 Govaere OR.csv")



head(genes)
sorted_pval <- genes[order(genes$pvalue),]
top_800 <- sorted_pval[1:800,]
dim(top_800)
write.csv(top_800, "top_800 Govaere OR.csv")



genes_normalized <- y$counts
top400_genes_normalized <- genes_normalized[rownames(top_400),]
top400_genes_normalized <- as.data.frame(top400_genes_normalized)

top400_genes_normalized
write.csv(top400_genes_normalized, "top400_genes_normalized for Gene ontology.csv")
save(top400_genes_normalized, file = 'top_400_genes')







