library("foreign")
library("np")
library("FNN")
library("corrplot")
library("pracma")
library("ggplot2")
library("glue")
#
vol_hypersphere <- function(r = 1, d = 1) {
    nom   = pi ^ (d/2)
    denom = gamma(d/2 + 1)
    unit_vol = nom / denom
    return ( unit_vol * (r ^ d) )
}

gauss_kernel <- function(x, sigma = 1) {
    # TODO: Should this assume that the kernel is centered at zero?

    const_term = 1 / (sqrt(2 * pi) * sigma)
    exp_term = exp( -(x ^ 2 / (2 * sigma ^ 2)))

    return  (const_term * exp_term)
}

BMP_truncated_inv_knn <- function(dataset, query, k = 2, t = 1) {
    # Breiman, Meisel, and Purcell - " From Location-Adaptive Density Estmiation Paper"
    # Notes: Volume and Radius is determined by each point in the dataset not the query.
    # Notes: Radius (used to calc Volume), is truncated .

    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)

    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(dataset,  dataset, k)[,k]
    radius_x = mapply(min, radius_x, (t * k) / n)
        # Addition: Adjust Estimates Where Duplication Causes Radius To Be Zero
        if (min(radius_x[radius_x > 0]) != Inf) {
            radius_x[radius_x == 0] = min(radius_x[radius_x > 0])
        }
        else  {
            radius_x[radius_x == 0] = 0.001
        }

    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)

    # Compute estimate at each query point
    out = c()
    for (p in query) {
        kern_out = gauss_kernel((p - dataset), radius_x)

        out = c(out, (1/n) * sum(kern_out))
    }
    return (out)
}
#
DATASET_NAME <- "Wilt Normalized No Dup."
FILE_PATH <- "datasets/semantic/Wilt/Wilt_withoutdupl_norm_05.arff"
N_DENSITY_POINTS <-  1024


df <-  read.arff(FILE_PATH)
df.inliers  <-  df[df$outlier == "no", ]
df.outliers <-  df[df$outlier == "yes", ]

# Each of the datasets should be
stopifnot("id" %in% names(df))
stopifnot("outlier" %in% names(df))

attribute_names <-  names(df)[ !names(df) %in% c("id", "outlier")]
for (attribute in attribute_names) {
    message(attribute)

    attribute_min <- min(df[, attribute])
    attribute_max <- max(df[, attribute])

    k.df      <- round(nrow(df) ^ (0.8))
    k.inlier  <- round(nrow(df.inliers) ^ (0.80))
    k.outlier <- round(nrow(df.outliers) ^ (0.80))

    density_query <- linspace(attribute_min, attribute_max, N_DENSITY_POINTS)

    density.overall  <- BMP_truncated_inv_knn(df[, attribute], density_query, k = k.df)
    density.inliers  <- BMP_truncated_inv_knn(df.inliers[, attribute], density_query, k = k.inlier)
    density.outliers <- BMP_truncated_inv_knn(df.outliers[, attribute], density_query, k = k.outlier)


    unique_values.inliers  <- length(unique(df.inliers[, attribute]))
    unique_values.outliers <- length(unique(df.outliers[, attribute]))

    ggplot() +
        geom_line(aes(x = density_query, y = density.overall, color = glue("All ({nrow(df)})"))) +
        geom_line(aes(x = density_query, y = density.inliers, color = glue("Inliers ({nrow(df.inliers)})"))) +
        geom_line(aes(x = density_query, y = density.outliers, color = glue("Outliers ({nrow(df.outliers)})"))) +
        geom_point(aes(x = df.inliers[, attribute], y = rep(0.01, nrow(df.inliers)), color = glue("Inliers ({nrow(df.inliers)})"))) +
        geom_point(aes(x = df.outliers[, attribute], y = rep(0.01, nrow(df.outliers)), color = glue("Outliers ({nrow(df.outliers)})"))) +
        scale_color_discrete(name = "Class") +
        ylab("Density Estimate") +
        xlab(glue("Normalized Attribute Range [{attribute_min},{attribute_max}]")) +
        ggtitle(label = glue("{DATASET_NAME}, {attribute} - Density Estimate of Inliers Vs Outliers")) +
        labs(subtitle = glue("At Zero: plot of values per class, unique inlier values = {unique_values.inliers}, unique outlier values = {unique_values.outliers}")) +
        labs(caption = glue("Overall K = {k.df}, Inlier K = {k.inlier}, Outlier K = {k.outlier}") ) # + scale_y_continuous(trans = "log10")

    ggsave(glue("Density-{DATASET_NAME}-{attribute}.png"))

}






