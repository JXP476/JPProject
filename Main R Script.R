install.packages("limma")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
install.packages("edgeR")
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
#NOT SURE IF GPLOTS ACTUALLY INSTALLED**
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("org.Dm.eg.db")
library(org.Dm.eg.db)



# FILTERING OUT LOWLY EXPRESSED GENES
  
counts <- geneExpressionCPM
targets <- patientData

table(patientData$Tissue)

mycpm <- cpm(geneExpressionCPM) # CY (22/12/2017): I SUSPECT THIS IS THE PROBLEM. YOU ARE APPLYING THE CPM FUNCTION TO THE COUNTS MATRIX. THIS FUNCTION DOES NOT APPLY HERE. THERE SHOULD BE ANOTHER MATRIX OF CPM VALUES (geneExpressionCPM?).
# Why isn't it on counts? don't understand why apply cpm function if already CPM?
head(mycpm)
thresh <- mycpm > 0.5
head(thresh)
table(rowSums(thresh))

keep <- rowSums(thresh) >= 7
table(keep)
counts.keep <- counts[keep,]
dim(counts.keep)



plot(counts[,1],mycpm[,1],xlim=c(0,15),ylim=c(0,5))
abline(v=10, col=2)


# CONVERT TO DGE LIST OBJECT 

y <- DGEList(counts.keep)
y
names(y)
y$samples
dim(y)
## ASK: Doesnt there need to be group 1 & group 2 (tumour vs normal)

# Library sizes and distribution plots

y$samples$lib.size
barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes")

logcpm <- cpm(y$counts,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# - COLOURED BY GROUPS:
group.col <- c("red","blue")[patientData$Tissue]
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)

# Multidimensional scaling plots

plotMDS(y)

# to make plot more informative; colour samples according to grouping info



# First we need a matrix of log counts:

logcounts <- cpm(y,log=TRUE)

# Hierarchical clustering with heatmap.2

logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
highly_variable_lcpm <- logcounts[select_var,]
dim(highly_variable_lcpm)
mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)

# Plot the heatmap- Normalisation
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=group.col,scale="row",margins=c(10,5))
y <- calcNormFactors(y)
y$samples
# THESE SHOW biased and unbiased MD plots side by side for the same sample to see the before and after TMM normalisation effect.
# P4T7
par(mfrow=c(1,2))
plotMD(logcounts,column=2)
abline(h=0,col="grey")
plotMD(y,column = 2)
abline(h=0,col="grey")

# P3T3
par(mfrow=c(1,2))
plotMD(logcounts,column=4)
abline(h=0,col="grey")
plotMD(y,column = 2)
abline(h=0,col="grey")

# P4N1
par(mfrow=c(1,2))
plotMD(logcounts,column=7)
abline(h=0,col="grey")
plotMD(y,column = 2)
abline(h=0,col="grey")

# Differential expression
#Set up design matrix
# We want to test for differences between the treated and untreated samples. However, we know that the library preparation adds variability to the data, so we need to account for it in our model. We do this by modelling both Group and Library as variables in our design matrix. This is known as an additive model.
design <- model.matrix(~patientData$Patient + patientData$Tissue)
design
colnames(design) <- c("Int","SEvsPE","UVsT")
par(mfrow=c(1,1))
v <- voom(y,design,plot=TRUE)
par(mfrow=c(1,2))
boxplot(logcounts)
abline(h=median(logcounts),col=4)
boxplot(v$E)
abline(h=median(v$E),col=4)