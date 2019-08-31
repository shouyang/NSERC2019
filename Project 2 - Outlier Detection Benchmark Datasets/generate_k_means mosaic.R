library(foreign)
library(glue)
library(cluster)
library(ggplot2)
library(dplyr)
library(reshape2)
library(treemap)
library(vcd)
#

DATASET_NAME <- "Wilt_withoutdupl_norm_05"
FILE_PATH <- "datasets/semantic/Wilt/Wilt_withoutdupl_norm_05.arff"
K_RANGE <- seq(2,15,1)






df <-  read.arff(FILE_PATH)
# df <- df[sample(1:48113, 30000),]


attribute_names   <-  names(df)[!names(df) %in% c("id", "outlier", "alphaLevel")]
attribute_indexes <-  subset(1:ncol(df), !names(df) %in% c("id", "outlier", "alphaLevel"))
#
cluster_cols    <- c()
silhouette_cols <- c()

for (i in K_RANGE) {
    temp_cluster_info <-  kmeans(df[, attribute_indexes], centers = i, nstart = 5, iter.max = 50)

    # Remap cluster results such that 1 is always the largest cluster by points and i is always the smallest by points
    sorted_clusters <- sort(table(temp_cluster_info$cluster), decreasing = T)
    sorted_clusters[1:i] <- 1:i


    new_col_name <- glue("K{i}")
    new_sil_name <- glue("Si{i}")
    df[, new_col_name]  <- recode(temp_cluster_info$cluster, !!!sorted_clusters)
    df[, new_sil_name] <-  silhouette(df[,new_col_name], dist(df[,attribute_names]))[, "sil_width"]

    cluster_cols    <- c(cluster_cols, new_col_name)
    silhouette_cols <- c(silhouette_cols, new_sil_name)
}

#
temp <- melt(df, id.vars = c("id", "outlier", attribute_names, silhouette_cols))

# temp2 <- table(temp$variable, temp$value, temp$outlier)
# mosaicplot(temp2, dir = "v", color = c("cadetblue1", "coral1", "red", "green"), xlab = "nstart = 100", ylab = "Cluster # (By Size)", las = 2, off = 10, main = DATASET_NAME)


struct <- structable(~ variable + value + outlier, data = temp, split_vertical = T)

labelling_table <- as.table(struct)
for (k_size in K_RANGE) {
    for (k_num in c(1,K_RANGE)) {
        cluster_at_k <- filter(temp, variable == glue("K{k_size}") & value == k_num)


        for (flag in c("yes", "no")) {
            sub_cluster_at_k <- filter(cluster_at_k, outlier == flag)
            if (nrow(sub_cluster_at_k) != 0) {
                sub_cluster_at_k_count      <- labelling_table[glue("K{k_size}"), k_num, flag]
                sub_cluster_at_k_mean_score <- round(mean(sub_cluster_at_k[, glue("Si{k_size}")]),2)

                # For the first cluster, append the global mean sil score.
                if (k_num == 1) {

                    k_mean_score_for_group = round(mean(filter(temp, variable == glue("K{k_size}") & outlier == flag)[,glue("Si{k_size}")]),2)

                    labelling_table[glue("K{k_size}"), k_num, flag] <- glue("{sub_cluster_at_k_count} | {sub_cluster_at_k_mean_score} \n  ASil:{k_mean_score_for_group}")
                }
                else {
                    labelling_table[glue("K{k_size}"), k_num, flag] <- glue("{sub_cluster_at_k_count} | {sub_cluster_at_k_mean_score}")
                }
            }
            else {
                if (k_num == 1) {
                    k_mean_score_for_group = round(mean(filter(temp, variable == glue("K{k_size}") & outlier == flag)[,glue("Si{k_size}")]),2)

                    labelling_table[glue("K{k_size}"), k_num, flag] <- glue("ASil:{k_mean_score_for_group}")
                }
                else {
                    labelling_table[glue("K{k_size}"), k_num, flag] <- " "
                }

            }
        }
    }
}


png( glue("{DATASET_NAME}-Kmeans-Mosaic-Silhouette-Plot.png"), width = 1920, height = 1080)
    mosaic(struct,
        main = glue("K-Means (2-15) Mosiac-Silhouette Plot of {DATASET_NAME}"),
        sub  = "Columns: Increasing K, Sorted By Cluster Size, Labels: Per Cluster # Obs | Sil Score, ASil = Mean Sil Score of class per k, Colour: Blue = Inliers, Red = Outliers",
        direction = "v",
        zero_size = 0,
        gp = gpar(fill = matrix( c("cadetblue1", "coral1"), 1, 2)),
        labeling = labeling_cells(text = labelling_table, gp_text = gpar(fontsize = 11, col = "black"), clip_cells = F))
dev.off()