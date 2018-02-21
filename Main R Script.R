install.packages("limma")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
install.packages("edgeR")
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
install.packages("gplots")
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("org.hs.eg.db")
source("https://bioconductor.org/biocLite.R")
biocLite("org.Hs.eg.db")
library(org.Hs.eg.db)

"load({Main R Script})"

rownames(geneExpressionCounts) = gsub(".*_(.*)", "\\1", row.names(geneExpressionCounts))
rownames(geneExpressionCPM) = gsub(".*_(.*)", "\\1", row.names(geneExpressionCPM))
rownames(geneExpressionTPM) = gsub(".*_(.*)", "\\1", row.names(geneExpressionTPM))


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
counts.keep <- geneExpressionCPM[keep,]
dim(counts.keep)

# Normalisation for composition bias

y <- DGEList(counts.keep)
y <- calcNormFactors(y)
y
y$samples
dim(y)


# Library sizes and distribution plots

y$samples$lib.size

barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes")

# logcpm

logcpm <- cpm(y$counts,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (normalised)")

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

# Adding annotation and saving the results - Add annotation from org.hs.eg.db (human samples)

columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)


# At the moment the row names are a concatenation of its HUGO Gene Symbol and Ensembl ID.
# THIS BIT OF CODE FOLLOWING TURNS ROW NAMES TO JUST Ensembl IDs
# This will replace the row names of the geneExpressionXXX matrix by the ensembl ID only:

rownames(geneExpressionCounts) = gsub(".*_(.*)", "\\1", row.names(geneExpressionCounts))
rownames(geneExpressionCPM) = gsub(".*_(.*)", "\\1", row.names(geneExpressionCPM))
rownames(geneExpressionTPM) = gsub(".*_(.*)", "\\1", row.names(geneExpressionTPM))

# This code then matches on the new row names:
# You can then filter based on Ensembl Gene IDs:

ann = select(org.Hs.eg.db, keys = rownames(fit),columns=c("ENSEMBL", "SYMBOL","GENENAME"),keytype="ENSEMBL")

head(ann)
ann

table(ann$ENSEMBL==rownames(fit))
fit$genes <- ann
topTable(fit,coef=3,sort.by="p")

# Alternative but didn't work:

ls("package:org.Hs.eg.db")

j <- toTable(org.Hs.egENSEMBL)
head(j)

symbol <- toTable(org.Hs.egSYMBOL)
genename <- toTable(org.Hs.egGENENAME)

ann1 <- merge(j,symbol,by="gene_id")
head(ann1)

ann2 <- merge(ann1,genename,by="gene_id")
head(ann2)

m <- match(rownames(fit),ann2$ensembl_id)
table(is.na(m)) 

ann3 <- ann2[m[!is.na(m)],]
head(ann3)
head(fit$genes)

topTable(fit,coef=3,sort.by="p")

# Checking expression of tetraspanin 6

ps <- grep("tetraspanin 6",fit$genes$GENENAME)
topTable(fit[ps,],coef=3)

out = sort(fit$p.value[, 3], decreasing=FALSE, index.return=TRUE)

fitSorted = fit
fitSorted$p.value = fit$p.value[out$ix, ]
fitSorted$genes = fit$genes[out$ix, ]

# PLOTS:

par(mfrow=c(1,2))

plotMD(fit,coef=3,status=results[,"Tumour/Normal"])

volcanoplot(fit,coef=3,highlight=100,names=fit$genes$SYMBOL)




vennDiagram # ***** need to figure out how to create this- good figure
