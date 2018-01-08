# ** 1) READING DATA INTO R **
install.packages("limma")
source("https://bioconductor.org/biocLite.R")
biocLite("limma")
library(limma)
install.packages("edgeR")
source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")
library(edgeR)
#NOT SURE IF GPLOTS ACTUALLY INSTALLED**
library(gplots)
install.packages("RColorBrewer")
library(RColorBrewer)
install.packages("org.Dm.eg.db")
library(org.Dm.eg.db)


# Initial quality control & data preparation:
# ** 2) FILTERING OUT LOWLY EXPRESSED GENES **
# (In test = TREATED VS UNTREATED. In mine = 1. TUMOUR VS NORMAL 2. T VS N in EACH PATIENT
# 1st interest is testing TUMOUR VS NORMAL groups
# To check how many samples we have each group= use the table command
counts <- geneExpressionCPM
targets <- patientData$Patient
# I think the above assignment step is irrelevant but have left in for now^
table(patientData$Tissue)
# ^ Normal-7, Tumour-20
# So the minimum sample size is 7
# Check the relationship between CPM and counts to see what CPM threshold we should be imposing
# Users should filter with CPMs rather than filtering on the counts directly
# latter does not account for differences in library sizes between samples
# Need to READ combine-australia general tutorial for info on how CPM filtering works
# Recall we’re looking for a CPM that corresponds to a count of roughly 10-15.
mycpm <- (geneExpressionCPM) # CY (22/12/2017): I SUSPECT THIS IS THE PROBLEM. YOU ARE APPLYING THE CPM FUNCTION TO THE COUNTS MATRIX. THIS FUNCTION DOES NOT APPLY HERE. THERE SHOULD BE ANOTHER MATRIX OF CPM VALUES (geneExpressionCPM?).
plot(counts[,1],mycpm[,1],xlim=c(0,20),ylim=c(0,50))
abline(v=10,col=2)
abline(h=2,col=4)
#^ looks at count data
# Which values in myCPM are greater than 0.5? (0.5 is CPM threshold)
thresh <- mycpm > 0.5
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 7
table(keep)
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- counts[keep,]
dim(counts.keep)
# = Filtering out 23869 out of total 39381?
# See whether threshold of 0.5 corresponds to a count of about 10-15 (it doesn't)
plot(mycpm[,1],counts[,1])
# Limit the x and y-axis to see what is happening at the smaller counts
plot(mycpm[,1],counts[,1],ylim=c(0,2),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)


# ** 3) CONVERT TO DGE LIST OBJECT ** 
# DGEList object holds the dataset to be analysed by edgeR and the subsequent calculations performed on the dataset
y <- DGEList(counts.keep) 

#^ Since the DGElist should contain : lib.size, norm.factors, group, genes


# ** 4) QUALITY CONTROL **
# Here conduct a number of quality QC plots. 1st, check the library sizes:
barplot(y$samples$lib.size)
# Next check the distribution of the counts using a boxplot:
par(mfrow=c(1,1))
# Get log2 counts per million- can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. The cpm function also adds a small offset to avoid taking log of zero.
logcpm <- cpm(y$counts,log=TRUE)
# Check distributions of samples using boxplots
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,outline=FALSE)
# Let's add a blue horizontal line that corresponds to the median logCPM
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# We can colour by our groups:
par(mfrow=c(1,2),oma=c(2,0,0,0))
group.col <- c("red","blue")[targets$Tissue]
# ^ shouldn't this colour the tumour and normal patients in red/blue accordingly?
boxplot(logcpm, xlab="", ylab="Log2 counts per million",las=2,col=group.col,
        pars=list(cex.lab=0.8,cex.axis=0.8))
abline(h=median(logcpm),col="blue")
title("Boxplots of logCPMs\n(coloured by groups)",cex.main=0.8)

# Any BIAS in data? P8N1 needs investigating further- looks completely diff

# Now check our MDS plots.
plotMDS(y)

# to make plot more informative; colour samples according to grouping info
plotMDS(y,col=col.cell)
n <-c("PT45", "PT44", "PT42")
Tissue <-c("Tumour", "Tumour", "Tumour")
col.cell <- c("purple","orange", "purple")
data.frame(sampleID,Tissue,col.cell)
legend("topleft",fill=c("purple","orange"),legend=levels(sampleID$Tissue))
levels(sampleID$Tissue)

#CURRENTLY CANT GET MDS TO WORK WITH COLOURED LABELS BUT ATTEMPTING TO FIX AGAINST TUTORIAL http://combine-australia.github.io/RNAseq-R/06-rnaseq-day1.html

# ** 5) Hierarchical clustering with heatmap **

# First we need a matrix of log counts:

logcounts <- cpm(y,log=TRUE)



