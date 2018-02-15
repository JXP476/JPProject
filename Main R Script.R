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
install.packages("org.hs.eg.db")
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

# Note: USE CAPITAL H in Hs


# FILTERING OUT LOWLY EXPRESSED GENES
  
counts <- geneExpressionCounts
targets <- patientData

table(patientData$Tissue)

mycpm <- geneExpressionCPM 
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

# logcpm

logcpm <- cpm(y$counts,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Multidimensional scaling plots

plotMDS(y)
col.cell <- c("purple","purple", "purple","purple","purple","purple","orange","orange","purple","purple","purple","purple","orange","purple","purple","purple","purple","orange","purple","purple","purple","purple","orange","purple","purple","orange","orange")
data.frame(patientData$Tissue,col.cell)
plotMDS(y,col=col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(patientData$Tissue))
title("Cell type")

# Hierarchical clustering with heatmap.2

logcounts <- cpm(y,log=TRUE)
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
highly_variable_lcpm <- logcounts[select_var,]
head(highly_variable_lcpm)
dim(highly_variable_lcpm)

mypalette <- brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
col.cell <- c("purple","purple", "purple","purple","purple","purple","orange","orange","purple","purple","purple","purple","orange","purple","purple","purple","purple","orange","purple","purple","purple","purple","orange","purple","purple","orange","orange")
data.frame(patientData$Tissue,col.cell)


# Plot the heatmap
heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")





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

# include the Patient + tissue as variables
# You are testing for differences between normal and tumour tissues whilst accounting for any systematic inter-patient differences

# SET UP DESIGN MATRIX
design <- model.matrix(~targets$Patient + targets$Tissue)
design
colnames(design) <- c("Intercept", "Patient","Tumour/Normal")

# Voom transform the data
par(mfrow=c(1,1))
v <- voom(y,design,plot=TRUE)
par(mfrow=c(1,2))
boxplot(logcounts)
abline(h=median(logcounts),col=4)
boxplot(v$E)
abline(h=median(v$E),col=4)


# Test for differential expression
fit <- lmFit(v,design)
fit <- eBayes(fit)
results <- decideTests(fit)
summary(results)

topTable(fit,coef=3,sort.by="p")

# Add annotation from org.hs.eg.db (human samples)

columns(org.Hs.eg.db)







vennDiagram # ***** need to figure out how to create this- good figure
