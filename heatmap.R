##############################
## R code to generate heatmap for RNAseq and RPPA data 
## using RSEM or HTseq counts for RNAseq data and DE
## RPPA data.
## Author: Sammed Mandape
## Contact: smandape@email.arizona.edu
## Date: Apr 16, 2019
##############################

# libraries
library(pheatmap)
library(dplyr)
library(gplots)
# Read data
data_input <- readRDS("6sample_rawCPM.RDS")
write.table(data_input, file= "6sample_rawCPM.csv", sep="\t")

# log2 normalization of data
logdata <- data_input
#logdata[,1:6] <- as.matrix(log2(logdata[1:6]))
logdataM <- as.matrix(log2(logdata[1:6]))

# checks
# head(logdataM)
# typeof(logdataM)
# dim(logdataM)
# sum(is.infinite(logdataM))
# sum(is.na(logdataM))

# # either remove Na/NaN/inf OR replace with 0
# to see if it has any INF value, can be also done for Na or NaN. Especially, when you get error something like 'Error in hclust(d, method = method) : 
# NA/NaN/Inf in foreign function call (arg 11)'. The sum should be zero.
### sum(is.infinite(data_subset_matrix))
#nrow(logdataM)
logdataM_sub <- logdataM
logdataM_sub[is.infinite(logdataM_sub)] <- 0
#head(logdataM_sub)

# remove rows that equal to 0 for all columns
# {x == 0 creates a matrix of logical values (TRUE/FALSE) and rowSums(x == 0) sums them up (TRUE == 1, FALSE == 0).
# Then you check if the sum of each row is not equal to the number of columns of your matrix (which are counted by ncol(x)).}
logdataM_sub <- logdataM_sub[rowSums(logdataM_sub == 0) != ncol(logdataM_sub),]
#nrow(logdataM_sub)
heatmap.2(logdataM_sub)
heatmap.2(logdataM_sub, scale="row", Rowv = FALSE, Colv = FALSE)
pheatmap(logdataM_sub)
write.table(logdataM_sub, file="logdataM_sub.txt", sep="\t")

