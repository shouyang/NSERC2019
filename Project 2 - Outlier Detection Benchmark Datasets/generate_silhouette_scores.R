library(foreign)
library(glue)
library(NbClust)
#
DATASET_NAME <- "WBC_withoutdupl_v10"
FILE_PATH <- "datasets/literature/WBC/WBC_withoutdupl_norm_v10.arff"


df <-  read.arff(FILE_PATH)
df.inliers  <-  df[df$outlier == "no", ]
df.outliers <-  df[df$outlier == "yes", ]

attribute_names   <-  names(df)[!names(df) %in% c("id", "outlier", "alphaLevel")]
attribute_indexes <-  subset(1:ncol(df), !names(df) %in% c("id", "outlier", "alphaLevel"))
#
NbClust(df[, attribute_names], max.nc = 100, method = "ward.D2", index = c("silhouette", "ch", "db"))
