# Initial setup and loading of required packages in R:

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
source("http://bioconductor.org/biocLite.R")
biocLite("GO.db")
library(GO.db)

"load({Main R Script})"
rownames(geneExpressionCounts) = gsub(".*_(.*)", "\\1", row.names(geneExpressionCounts))
rownames(geneExpressionCPM) = gsub(".*_(.*)", "\\1", row.names(geneExpressionCPM))
rownames(geneExpressionTPM) = gsub(".*_(.*)", "\\1", row.names(geneExpressionTPM))


# 1) Filtering of lowly expressed genes
  
counts <- geneExpressionCounts
targets <- patientData

table(patientData$Tissue)

mycpm <- geneExpressionCPM 
head(mycpm)
thresh <- mycpm > 0.5
head(thresh)
table(rowSums(thresh))

keep <- rowSums(thresh) >= 2
table(keep)
counts.keep <- geneExpressionCPM[keep,]
dim(counts.keep)

# 2) TMM Normalisation

y <- DGEList(counts.keep)
y <- calcNormFactors(y)
y$samples
y

# P5N2 - Before and after TMM Normalisation
par(mfrow=c(1,2))
plotMD(logcounts,column=13)
abline(h=0,col="red")

plotMD(y,column = 13)
abline(h=0,col="red")

# Effective library sizes and distribution plots

y$samples$lib.size

barplot(y$samples$lib.size,names=colnames(y),las=2)
title("Barplot of library sizes (Normalised)")

barplot(patientData$total_counts,names=colnames(y),las=2)
title("Barplot of library sizes (Unnormalised)")

# logcpm

logcpm <- cpm(y$counts,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (Normalised)")

# Dr Yau- please could you let me know if this following 77-80 is correct for an unnormalised boxplot?
logcpm <- cpm(counts,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (Unnormalised)")


# Multidimensional scaling plots

plotMDS(y)
col.cell <- c("purple","purple", "purple","purple","purple","purple","orange","orange","purple","purple","purple","purple","orange","purple","purple","purple","purple","orange","purple","purple","purple","purple","orange","purple","purple","orange","orange")
data.frame(patientData$Tissue,col.cell)
plotMDS(y,col=col.cell)
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

heatmap.2(highly_variable_lcpm,col=rev(morecols(50)),trace="none", main="Top 500 most variable genes across samples",ColSideColors=col.cell,scale="row")

# Differential expression

# Set up design matrix
design <- model.matrix(~targets$Patient + targets$Tissue)
design
colnames(design) <- c("Intercept", "Patient","Tumour/Normal")

# Voom transform the data
par(mfrow=c(1,1))
v <- voom(y,design,plot=TRUE)
v

# Test for differential expression
fit <- lmFit(v,design)
fit <- eBayes(fit)
results <- decideTests(fit)
results
summary(results)
topTable(fit,coef=3,sort.by="p") # **** Does this need to be here?

# Adding annotation and saving the results

columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)

ann = select(org.Hs.eg.db, keys = rownames(fit),columns=c("ENSEMBL", "SYMBOL","GENENAME"),keytype="ENSEMBL")
fit$genes <- ann

tmp = topTable(fit,coef="Tumour/Normal", sort.by="p", number=30)
cat(paste(tmp$SYMBOL, collapse = '\n'))

# PLOTS after testing for DE

par(mfrow=c(1,2))

plotMD(fit,coef=3,status=results[,"Tumour/Normal"])

volcanoplot(fit,coef=3,highlight=100,names=fit$genes$SYMBOL)

tmp = topTable(fit,coef="Tumour/Normal", sort.by="p", number=30)
cat(paste(tmp$SYMBOL, collapse = '\n'))
tmp$SYMBOL

# Gene set testing with Goana

go <- goana(tmp$SYMBOL)
topGO(go, n=15)
head(tmp)

keg <- kegga(tmp$SYMBOL, species="Hs")


topKEGG(keg, n=15, truncate=34)


