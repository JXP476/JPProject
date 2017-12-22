library(limma)
library(edgeR)
library(gplots)
library(RColorBrewer)
library(org.Dm.eg.db)
# Initial quality control & data preparation:
# ** 1st= Filtering out lowly expressed genes **
# (In test = TREATED VS UNTREATED. In mine = 1. TUMOUR VS NORMAL 2. T VS N in EACH PATIENT
# 1st interest is testing TUMOUR VS NORMAL groups
# To check how many samples we have each group= use the table command
counts <- geneExpressionCPM
targets <- patientData
# I think the above assignment step is irrelevant but have left in for now^
table(patientData$Tissue)
# ^ Normal-7, Tumour-20
# So the minimum sample size is 7
# Check the relationship between CPM and counts to see what CPM threshold we should be imposing
# Users should filter with CPMs rather than filtering on the counts directly
# latter does not account for differences in library sizes between samples
# Need to READ combine-australia general tutorial for info on how CPM filtering works
# Recall we’re looking for a CPM that corresponds to a count of roughly 10-15.
mycpm <- cpm(counts)
plot(counts[,1],mycpm[,1],xlim=c(0,20),ylim=c(0,50))
abline(v=10,col=2)
abline(h=2,col=4)
#^ looks at count data
# Which values in myCPM are greater than 0.5? (0.5 is CPM threshold)
thresh <- mycpm > 0.5
table(rowSums(thresh))
# we would like to keep genes that have at least 2 TRUES in each row of thresh
keep <- rowSums(thresh) >= 2
table(keep)
# Subset the rows of countdata to keep the more highly expressed genes
counts.keep <- counts[keep,]
summary(keep)
dim(counts.keep)
# Filter out 15512 ?
# See whether threshold of 0.5 corresponds to a count of about 10-15
plot(mycpm[,1],counts[,1])
# Limit the x and y-axis to see what is happening at the smaller counts
plot(mycpm[,1],counts[,1],ylim=c(0,2),xlim=c(0,3))
# Add a vertical line at 0.5 CPM
abline(v=0.5)
abline(v=0.5)
j

