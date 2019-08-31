library(foreign)
library(GGally)
library(glue)
#
DATASET_NAME <- "Wilt_withoutdupl_norm_05"
FILE_PATH <- "datasets/semantic/Wilt/Wilt_withoutdupl_norm_05.arff"


df <-  read.arff(FILE_PATH)
df.inliers  <-  df[df$outlier == "no", ]
df.outliers <-  df[df$outlier == "yes", ]


df$alphaLevel[df$outlier == "yes"] <- 0.4
df$alphaLevel[df$outlier == "no" ] <- 0.1

attribute_indexes <-  subset(1:ncol(df), !names(df) %in% c("id", "outlier", "alphaLevel"))


ggparcoord(df, columns = attribute_indexes, groupColumn = "outlier", alphaLines = "alphaLevel") +
    ggtitle(label = glue("{DATASET_NAME}, {attribute} - Density Estimate of Inliers Vs Outliers")) +
    labs(subtitle = glue("Inliers = {nrow(df.inliers)}, Outliers = {nrow(df.outliers)}"))

ggsave(glue("{DATASET_NAME} Parallel Plot.png"))
