library(foreign)
library(glue)
library(cluster)
library(ggplot2)
library(dplyr)
library(reshape2)
library(treemap)
#

DATASET_NAME <- "WPBC_withoutdupl_norm"
FILE_PATH <- "datasets/semantic/WPBC/WPBC_withoutdupl_norm.arff"
K_RANGE <- seq(2,15,1)






df <-  read.arff(FILE_PATH)
df.inliers  <-  df[df$outlier == "no", ]
df.outliers <-  df[df$outlier == "yes", ]

attribute_names   <-  names(df)[!names(df) %in% c("id", "outlier", "alphaLevel")]
attribute_indexes <-  subset(1:ncol(df), !names(df) %in% c("id", "outlier", "alphaLevel"))
#


for (i in K_RANGE) {
    temp_cluster_info <-  kmeans(df[, attribute_indexes], centers = i, nstart = 10, iter.max = 100)

    # Remap cluster results such that 1 is always the largest cluster by points and i is always the smallest by points
    sorted_clusters <- sort(table(temp_cluster_info$cluster), decreasing = T)
    sorted_clusters[1:i] <- 1:i

    df[, glue("K = {i}")] <- recode(temp_cluster_info$cluster, !!!sorted_clusters)
}

temp <- melt(df, id.vars = c("id", "outlier", attribute_names))
temp2 <- table(temp$variable, temp$value, temp$outlier)


png(glue("{DATASET_NAME} - mosaic plot.png"), width=1366, height=768)
    mosaicplot(temp2, dir = "v", color = c("cadetblue1", "coral1"), xlab = "nstart = 100", ylab = "Cluster # (By Size)", las = 2, off = 10, main = DATASET_NAME)
dev.off()
