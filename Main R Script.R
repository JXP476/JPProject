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
  
counts <- geneExpressionCounts
targets <- patientData$Patient

table(patientData$Tissue)

mycpm <- (geneExpressionCPM) # CY (22/12/2017): I SUSPECT THIS IS THE PROBLEM. YOU ARE APPLYING THE CPM FUNCTION TO THE COUNTS MATRIX. THIS FUNCTION DOES NOT APPLY HERE. THERE SHOULD BE ANOTHER MATRIX OF CPM VALUES (geneExpressionCPM?).
head(mycpm)
thresh <- mycpm > 1
head(thresh)
table(rowSums(thresh))

keep <- rowSums(thresh) >= 7
table(keep)
counts.keep <- counts[keep,]
dim(counts.keep)

plot(mycpm[,1],counts[,1])

# CONVERT TO DGE LIST OBJECT 
# DGEList object holds the dataset to be analysed by edgeR and the subsequent calculations performed on the dataset
# Since the DGElist should contain : lib.size, norm.factors, group, genes
y <- DGEList(counts.keep) 
y
names(y)
y$samples


# ** 4) QUALITY CONTROL **
# Now that we have got rid of the lowly expressed genes and have our counts stored in a DGEList object
# Look at plots to check that the data is good quality, and that the samples are as we would expect.
# Here conduct a number of quality QC plots

# 1st- Library sizes and distribution plots
y$samples$lib.size
barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes")
# Count data is not normally distributed, so if we want to examine the distributions of the raw counts we need to log the counts

# Next check the distribution of the counts using a boxplot:
# Get log2 counts per million- can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. The cpm function also adds a small offset to avoid taking log of zero.
logcpm <- cpm(y$counts,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# - COLOURED BY GROUPS:
group.col <- c("red","blue")[patientData$Tissue]
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)

# Any BIAS in data? P8N1 needs investigating further- looks completely diff

# 2nd- Multidimensional scaling plots

plotMDS(y)

# to make plot more informative; colour samples according to grouping info
plotMDS(y,col=col.cell)
data.frame(sampleID,Tissue,col.cell)
sampleID <-c("P4T5", "P4T7", "P3T1", "P3T3", "P3T5", "P3T7", "P4N1", "P4N2", "P4T1", "P4T3", "P1T2", "P1T6", "P5N2", "P5T2", "" )

#CURRENTLY CANT GET MDS TO WORK WITH COLOURED LABELS BUT ATTEMPTING TO FIX AGAINST TUTORIAL http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html

# ** 5) Hierarchical clustering with heatmap **
sam
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