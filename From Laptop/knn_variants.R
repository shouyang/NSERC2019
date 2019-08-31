library("FNN")


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


BMP_knn <- function(dataset, query, k = 2) {
    # Breiman, Meisel, and Purcell - " From Location-Adaptive Density Estmiation Paper"
    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)
    
    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(test_data,  test_data, k)[,k]
    
    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)

    # Compute estimate at each query point
    out = c()
    for (p in query) {
        kern_out = gauss_kernel((p - test_data), radius_x)
        
        out = c(out, (1/n) * sum(kern_out))
    }
    
    return (out)
}

BMP_truncated_knn <- function(dataset, query, k = 2, t = 1) {
    # Breiman, Meisel, and Purcell - " From Location-Adaptive Density Estmiation Paper"
    n = dim(t(t(dataset)))[1]
    d = dim(t(t(dataset)))[2]
    k = round(k)
    
    # Get each data point's k-th nearest neighbour distance (radius)
    radius_x = knnx.dist(test_data,  test_data, k)[,k]
    radius_x = min(radius_x, (t * k) / n)
    
    # Compute the two compoents dependent on each data point
    volume_x = vol_hypersphere(radius_x, d)
    
    # Compute estimate at each query point
    out = c()
    for (p in query) {
        kern_out = gauss_kernel((p - test_data), radius_x)
        
        out = c(out, (1/n) * sum(kern_out))
    }
    
    return (out)
}



library(benchden)


N_SAMPLES  = 1000
test_data  = rberdev(N_SAMPLES, 24)
EVAL_RANGE = seq(-3,3,0.01)
k  = 10


BMP_estimate       = BMP_knn(test_data, EVAL_RANGE, k = N_SAMPLES ^ (4/6))
truncated_estmiate = BMP_truncated_knn(test_data, EVAL_RANGE, k = N_SAMPLES ^ (4/6))

plot(EVAL_RANGE, dberdev(EVAL_RANGE, 24), type = "l")
    lines(EVAL_RANGE, BMP_estimate, col = "red")
    lines(EVAL_RANGE, truncated_estmiate, col = "blue") 