library(foreign)
library(glue)
library(dplyr)
library(calibrateBinary)
library(caret)
# Import Dataset
DATASET_NAME <- "Shuttle_withoutdupl_norm_v10"
FILE_PATH <- "datasets/literature/Shuttle/Shuttle_withoutdupl_norm_v10.arff"
K_RANGE <- seq(2,15,1)

df <-  read.arff(FILE_PATH)
df$outlier <- recode(df$outlier, yes = 0, no = 1)

attribute_names   <-  names(df)[!names(df) %in% c("id", "outlier")]
attribute_indexes <-  subset(1:ncol(df), !names(df) %in% c("id", "outlier"))


dem <- length(attribute_names)




df.train <- data.matrix(df[c(4:1000),])
df.test  <- data.matrix(df[1:3,])

cv_res <- cv.KLR(df.train[,attribute_indexes], df.train[,"outlier"])
scores <- KLR(df.train[,attribute_indexes], df.train[,"outlier"], df.test[, attribute_indexes], lambda = cv_res$lambda, rho = cv_res$rho)



