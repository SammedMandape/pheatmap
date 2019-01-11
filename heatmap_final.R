######### 
## Author: Sammed N. Mandape 
## This R code is to generate heatmap for expression values, specifically this code plots a heatmap for FPKM values.
## Scaling was done using zFPKM library. 
## FYI: Codes are commented by '###', if you want to run this line. 
#########

## Read data
data_input <- read.delim('Genes.csv', header = T, row.names = "Genes", sep=",")

# scaling data specific for FPKM values using zFPKM library, install zFPKM library if needed.
### BiocManager::install("zFPKM", version = "3.8")
library('zFPKM')
data_subset <- zFPKM(data_input[,1:4])

# convert to matrix
data_subset_matrix <- as.matrix(data_subset)

# either remove Na/NaN/inf OR replace with 0
# to see if it has any INF value, can be also done for Na or NaN. Especially, when you get error something like 'Error in hclust(d, method = method) : 
# NA/NaN/Inf in foreign function call (arg 11)'. The sum should be zero.
### sum(is.infinite(data_subset_matrix))
data_subset_matrix[is.infinite(data_subset_matrix)] <- 0

# plot heatmap, here pheatmap as well as heatmap.2 is used for comaprison
pheatmap(data_subset_matrix)
heatmap.2(data_subset_matrix)

# the following is for adding row and column annotations to the heatmap
my_hclust_gene<-hclust(dist(data_subset_matrix), method="complete")
# install if necessary
install.packages("dendextend")
# load package
library(dendextend)

as.dendrogram(my_hclust_gene) %>%
  plot(horiz = TRUE)

# define number of clusters you want in 'k'
my_gene_col <- dendextend::cutree(tree = as.dendrogram(my_hclust_gene), k = 6)

# format cluster names into a separate column and put it in a dataframe
my_gene_col1 <- data.frame(cluster = ifelse(my_gene_col == 1, "Cluster 1", ifelse(my_gene_col == 2, "Cluster 2", ifelse(my_gene_col == 3, "Cluster 3", ifelse(my_gene_col == 4, "Cluster 4", ifelse(my_gene_col== 5, "Cluster 5", "Cluster 6"))))))

# adding a second row annotation from data_input, the second annotation is genes grouped into functional classes
my_gene_col2 <- merge(my_gene_col1, data_input, by=0)

# define rownames 
rownames(my_gene_col2)=my_gene_col2$Row.names

# rename column names
my_gene_anno<-my_gene_col2[,c("cluster","Groups")]

# plot heatmap and change font size
pheatmap(data_subset_matrix,annotation_row = my_gene_anno)
pheatmap(data_subset_matrix,annotation_row = my_gene_anno, fontsize_row = 8, fontsize_col = 10)

#customizing colors of cluster legend
Clusters = c("#1B9E77","firebrick","#7570B3","#E7298A","#66A61E","#00008B")
names(Clusters) = c("Cluster 1","Cluster 2", "Cluster 3", "Cluster 4","Cluster 5","Cluster 6")
ann_colors = list(Clusters = Clusters)
pheatmap(data_subset_matrix, fontsize_row = 8, fontsize_col = 10, annotation_row = my_gene_anno, annotation_colors = ann_colors)


# here,grouping / ordering is done just by functioanl classes, second row annotation defined above
# 1. heatmap is with clustering but wont' show the actual tree (treeheight_row = F)
# 2. heatmap is without clustering or it will be grouped together by functional classes defined.
pheatmap(data_subset_matrix,annotation_row = my_gene_anno)
rownames(my_gene_anno1)=rownames(data_input)
pheatmap(data_subset_matrix,annotation_row = my_gene_anno1, clustering_distance_rows = "euclidean",clustering_method = "complete",treeheight_row = F)
pheatmap(data_subset_matrix,annotation_row = my_gene_anno1, cluster_row = F)


##ignore the following

# matrix withour zFPKM normalization or scaling
data_input_subset <- as.matrix(data_input[,1:4])
# heatmap with row scaling
heatmap.2(data_input_subset, scale = "row")