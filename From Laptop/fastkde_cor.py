#!python

import numpy as np
import scipy.stats as scipy
from fastkde import fastKDE
import pylab as PP

# Generate two random variables dataset (representing 100000 pairs of datapoints)
N_options = [100,1000,10000]
RUNS_PER_N = 100
cor = {}

for N in N_options:
    cor[N] = []
    for run in range(RUNS_PER_N):
        test_data  = np.random.lognormal(size = N)

        # Do the self-consistent density estimate

        x = [np.linspace(-5,5, 513)]
        myPDF, x_axis = fastKDE.pdf(test_data, axes = x )

        # Construct filter for low density regions 
        low_density = scipy.lognorm(1).cdf(x_axis)

        lower_low_density =  scipy.lognorm(1).cdf(x_axis) < 1
        upper_low_density =  (1- scipy.lognorm(1).cdf(x_axis)) < 1

        low_density = [lower_low_density | upper_low_density]

        cor[N].append( scipy.spearmanr(myPDF[low_density], scipy.lognorm(1).pdf(x_axis)[low_density])[0] )

    PP.plot(x_axis, myPDF, color = "red", label = "estimate" )
    PP.plot(test_data, [0.01] * len(test_data), "+", color = "red", label = "test_data n = {0}".format(N))
    PP.plot(x_axis, scipy.lognorm(1).pdf(x_axis), color = "black", label = "actual")
    PP.legend(loc= "upper left")
    PP.show()

PP.boxplot( (cor[100],cor[1000],cor[10000]) ) 
PP.xticks([1,2,3],["100","1000","10000"])
PP.title("Low density corr of fastkde on Normal Distribution - Various n")
PP.show()
