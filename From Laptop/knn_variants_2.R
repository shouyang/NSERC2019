library("FNN")

vol_hypersphere <- function(r = 1, d = 1) {
    nom   = pi ^ (d/2)
    denom = gamma(d/2 + 1)
    unit_vol = nom / denom
    return ( unit_vol * (r ^ d) )
}

gauss_kernel <- function(x, sigma = 1, d = 1) {
    # TODO: Should this assume that the kernel is centered at zero?
    
    const_term =  1 / (sqrt(2 * pi) * sigma) ^ d 
    exp_term = exp((-1 / 2) * ( norm(x, type = "2") ^ 2 / sigma ^ 2) )
    
    return  (const_term * exp_term)
}

kernel_knn <- function(dataset, query, k = 2) {
    # Basic KNN density estimation as defined by - " From Location-Adaptive Density Estmiation Paper"
    # Notes: Volume and Radius of data point per query is fixed by k-th nearest neighbours.
    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)
    
    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(test_data,  query, k)[,k]
    
    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)
    
    # Compute estimate at each query point
    out = c()
    for (i in 1:nrow(query)) {
        row_query = matrix(query[i,], nrow = 1)
        
        kern_out =  apply((row_query - test_data), 1, gauss_kernel, sigma = radius_x[i], d = d)
        out = c(out, (1 / n) * sum(kern_out))
    }
    
    return (out)
}

BMP_knn <- function(dataset, query, k = 2) {
    # Breiman, Meisel, and Purcell- " From Location-Adaptive Density Estmiation Paper"
    # Notes: Volume and Radius is determined by each point in the dataset not the query. 
    
    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)
    
    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(test_data,  test_data, k)[,k]
    
    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)
    
    # Compute estimate at each query point
    out = c()
    for (i in 1:nrow(query)) {
        row_query = matrix(query[i,], nrow = 1)
        
        kern_out = c()
        for (j in 1:nrow(dataset)) {
            temp = row_query - dataset[j,]
            temp = gauss_kernel(temp, radius_x[j], d)

            kern_out = c(kern_out, temp)
        }
        
        out = c(out, (1/n) * sum(kern_out))
    }
    return (out)
}

BMP_truncated_knn <- function(dataset, query, k = 2, t = 1) {
    # Breiman, Meisel, and Purcell - " From Location-Adaptive Density Estmiation Paper"
    # Notes: Volume and Radius is determined by each point in the dataset not the query. 
    # Notes: Radius (used to calc Volume), is truncated .
    
    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)
    
    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(test_data,  test_data, k)[,k]
    radius_x = mapply(min, radius_x, (t * k) / n)
    
    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)
    
    # Compute estimate at each query point
    out = c()
    for (i in 1:nrow(query)) {
        row_query = matrix(query[i,], nrow = 1)
        
        kern_out = c()
        for (j in 1:nrow(dataset)) {
            temp = row_query - dataset[j,]
            temp = gauss_kernel(temp, radius_x[j], d)
            
            kern_out = c(kern_out, temp)
        }
        
        out = c(out, (1/n) * sum(kern_out))
    }
    return (out)
}

BMP_truncated_inv_knn <- function(dataset, query, k = 2, t = 4) {
    # Breiman, Meisel, and Purcell - " From Location-Adaptive Density Estmiation Paper"
    # Notes: Volume and Radius is determined by each point in the dataset not the query. 
    # Notes: Radius (used to calc Volume), is truncated .
    
    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)
    
    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(test_data,  test_data, k)[,k]
    radius_x = mapply(max, radius_x, (t * k) / n)
    
    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)
    
    # Compute estimate at each query point
    out = c()
    for (i in 1:nrow(query)) {
        row_query = matrix(query[i,], nrow = 1)
        
        kern_out = c()
        for (j in 1:nrow(dataset)) {
            temp = row_query - dataset[j,]
            temp = gauss_kernel(temp, radius_x[j], d)
            
            kern_out = c(kern_out, temp)
        }
        
        out = c(out, (1/n) * sum(kern_out))
    }
    return (out)
}

library(ggplot2)
library(benchden)
library(np)


N_SAMPLES = 100
K = round(N_SAMPLES ^ (4/5))
t = 1
    
    
test_data = data.frame( a = rnorm(N_SAMPLES), b = rnorm(N_SAMPLES))

x_query = seq(-5,5, 0.5)
y_query = seq(-5,5, 0.5)
query = expand.grid(x = x_query, y = y_query)

kernel_knn_estmiate             = kernel_knn(test_data, query, k = K)
BMP_knn_estmiate                = BMP_knn(test_data, query, k = K)
BMP_truncated_knn_estmiate      = BMP_truncated_knn(test_data, query, k = K, t)
BMP_truncated_inv_knn_estimate  = BMP_truncated_inv_knn(test_data, query, k = K, t)


plot1 <- ggplot(cbind(query, kernel_knn_estmiate), aes(x = x, y = y, z = kernel_knn_estmiate)) +
    stat_contour(geom = "polygon", aes(fill =..level..)) +
    geom_tile(aes(fill= kernel_knn_estmiate)) +
    stat_contour(bins = 30) +
    ggtitle("Kernel Knn Estmiate") + 
    labs( subtitle = sprintf("n = %d, k = %d", N_SAMPLES, K))

print(plot1)

plot2 <- ggplot(cbind(query, BMP_knn_estmiate), aes(x = x, y = y, z = BMP_knn_estmiate)) +
    stat_contour(geom = "polygon", aes(fill =..level..)) +
    geom_tile(aes(fill= BMP_knn_estmiate)) +
    stat_contour(bins = 30) +
    ggtitle("BMP Knn Estmiate") + 
    labs( subtitle = sprintf("n = %d, k = %d", N_SAMPLES, K))

print(plot2)

plot3 <- ggplot(cbind(query, BMP_truncated_knn_estmiate), aes(x = x, y = y, z = BMP_truncated_knn_estmiate)) +
    stat_contour(geom = "polygon", aes(fill =..level..)) +
    geom_tile(aes(fill= BMP_truncated_knn_estmiate)) +
    stat_contour(bins = 30) +
    ggtitle("BMP Knn w/ Truncation Estmiate") +
    labs( subtitle = sprintf("n = %d, k = %d, t = %d", N_SAMPLES, K, t))


print(plot3)

plot4 <- ggplot(cbind(query, BMP_truncated_inv_knn_estimate), aes(x = x, y = y, z = BMP_truncated_inv_knn_estimate)) +
    stat_contour(geom = "polygon", aes(fill =..level..)) +
    geom_tile(aes(fill= BMP_truncated_inv_knn_estimate)) +
    stat_contour(bins = 30) +
    stat_contour(bins = 30, data = cbind(query, BMP_knn_estmiate), aes(x = x, y = y, z = BMP_knn_estmiate, color = "BMP Knn")) +
    stat_contour(bins = 30, data = cbind(query, BMP_truncated_knn_estmiate), aes(x = x, y = y, z = BMP_truncated_knn_estmiate, color = "BMP Truncated Knn")) +
    
    ggtitle("BMP Knn w/ Inv. Truncation Estmiate") +
    labs( subtitle = sprintf("n = %d, k = %d, t = %d", N_SAMPLES, K, t))


print(plot4)






