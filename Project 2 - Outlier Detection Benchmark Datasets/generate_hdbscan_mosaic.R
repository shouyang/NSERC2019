library(foreign)
library(glue)
library(cluster)
library(ggplot2)
library(dplyr)
library(reshape2)
library(dbscan)
library(vcd)
#

DATASET_NAME <- "Ionosphere_withoutdupl_norm"
FILE_PATH <- "datasets/literature/Ionosphere/Ionosphere_withoutdupl_norm.arff"
K_RANGE <- seq(2,15,1)

df <-  read.arff(FILE_PATH)
# df <- df[sample(1:48113, 30000),]

attribute_names   <-  names(df)[!names(df) %in% c("id", "outlier", "alphaLevel")]
attribute_indexes <-  subset(1:ncol(df), !names(df) %in% c("id", "outlier", "alphaLevel"))
#


min_pts <- 7
res <-  hdbscan(df[, attribute_names], min_pts, gen_hdbscan_tree = T, gen_simplified_tree = F)

# Extract meaningful eps values from hdbscan tree
useful_eps = rep(-1,15)
cluster_cols    <- c()
silhouette_cols <- c()


cur_eps = 1

candidate_eps <- sort(res$hc$height, decreasing = T)

if (length(candidate_eps) > 100) {
    candidate_eps <- candidate_eps[ c(1:100, seq(100, length(candidate_eps), 5))]
}


for (eps in candidate_eps)
{
    print( c(cur_eps, length(candidate_eps)))
    cur_eps <- cur_eps +  1

    x <-  dbscan(df[, attribute_names], eps - 0.0001, min_pts, borderPoints = F)
    k <- length(unique(x$cluster[x$cluster != 0]))


    if ( k >= 2 && k <= 15 && useful_eps[k] == -1) {
        useful_eps[k] <- eps - 0.0001

        new_col_name <- glue("K{k}")
        new_sil_name <- glue("Si{k}")

        sorted_clusters <- sort(table(x$cluster[x$cluster != 0]), decreasing = T)
        sorted_clusters[1:k] <- 1:k
        sorted_clusters["0"] <- 0

        df[new_col_name] <- recode(x$cluster, !!!sorted_clusters)
        df[new_sil_name] <- silhouette(df[, new_col_name], dist(df[, attribute_names]))[, "sil_width"]

        max_k <- k
        cluster_cols    <- c(cluster_cols, new_col_name)
        silhouette_cols <- c(silhouette_cols, new_sil_name)
    }

    #  Stop when 15 clusters are found for large list of eps values.
    if (useful_eps[15] != -1) {
        break
    }
}

#
#
temp <- melt(df, id.vars = c("id", "outlier", attribute_names, silhouette_cols))

# temp2 <- table(temp$variable, temp$value, temp$outlier)
# mosaicplot(temp2, dir = "v", color = c("cadetblue1", "coral1", "red", "green"), xlab = "nstart = 100", ylab = "Cluster # (By Size)", las = 2, off = 10, main = DATASET_NAME)


struct <- structable(~ variable + value + outlier, data = temp, split_vertical = T)

labelling_table <- as.table(struct)
for (k_size in 2:max_k) {
    for (k_num in 0:max_k) {
        cluster_at_k <- filter(temp, variable == glue("K{k_size}") & value == k_num)

        if (nrow(cluster_at_k) == 0) {
            next
        }


        for (flag in c("yes", "no")) {
            sub_cluster_at_k <- filter(cluster_at_k, outlier == flag)
            if (nrow(sub_cluster_at_k) != 0) {
                sub_cluster_at_k_count      <- labelling_table[glue("K{k_size}"), k_num + 1, flag]
                sub_cluster_at_k_mean_score <- round(mean(sub_cluster_at_k[, glue("Si{k_size}")]),2)

                # For the first cluster, append the global mean sil score.
                if (k_num == 0) {

                    k_mean_score_for_group = round(mean(filter(temp, variable == glue("K{k_size}") & outlier == flag)[,glue("Si{k_size}")]),2)

                    if (flag == "no") {
                        labelling_table[glue("K{k_size}"), k_num + 1, flag] <- glue("(Noise Cluster) \n {sub_cluster_at_k_count} | {sub_cluster_at_k_mean_score} \n ASil: {k_mean_score_for_group} \n EPS: {round(useful_eps[k_size],4)}")
                    }
                    else {
                        labelling_table[glue("K{k_size}"), k_num + 1, flag] <- glue("{sub_cluster_at_k_count} | {sub_cluster_at_k_mean_score} \n ASil: {k_mean_score_for_group}")
                    }
                }
                else {
                    labelling_table[glue("K{k_size}"), k_num + 1, flag] <- glue("{sub_cluster_at_k_count} | {sub_cluster_at_k_mean_score}")
                }
            }
            else {
                if (k_num == 0) {
                    k_mean_score_for_group = round(mean(filter(temp, variable == glue("K{k_size}") & outlier == flag)[,glue("Si{k_size}")]),2)

                    labelling_table[glue("K{k_size}"), k_num + 1, flag] <- glue("ASil: {k_mean_score_for_group}")
                }

            }
        }
    }
}


png( glue("{DATASET_NAME}-HDBSCAN-Mosaic-Silhouette-HB-Plot.png"), width = 1920, height = 1080)
mosaic(struct,
       main = glue("HDBSCAN (K = 2-N, Max = 15) Mosiac-Silhouette Plot of {DATASET_NAME}"),
       sub  = "Columns: Increasing K, Sorted By Cluster Size (Except Noise Cluster First), Labels: Per Cluster # Obs | Sil Score, ASil = Mean Sil Score of class per k EPS = Cut Threshold EPS, Colour: Blue = Inliers, Red = Outliers",
       direction = "v",
       zero_size = 0,
       gp = gpar(fill = matrix( c("cadetblue1", "coral1"), 1, 2)),
       labeling = labeling_cells(text = labelling_table, gp_text = gpar(fontsize = 11, col = "black"), clip_cells = F))
dev.off()

