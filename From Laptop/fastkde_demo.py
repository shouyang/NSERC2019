#!python

import numpy as np
from fastkde import fastKDE
import pylab as PP

#Generate two random variables dataset (representing 100000 pairs of datapoints)
N = 1000
mu, sigma = 0., 1. # mean and standard deviation
var1 = np.random.lognormal(mu,sigma, size=int(N))

#Do the self-consistent density estimate
myPDF,axes = fastKDE.pdf(var1)
#Plot contours of the PDF should be a set of concentric ellipsoids centered on
#(0.1, -300) Comparitively, the y axis range should be tiny and the x axis range
#should be large


pdf = (np.exp(-(np.log(axes) - mu)**2 / (2 * sigma**2)) / (axes * sigma * np.sqrt(2 * np.pi)))


print(min(var1))


PP.plot(axes,myPDF)
PP.plot(axes, pdf, color="r")
PP.plot(var1, 1000 * [0.1], "ro")

PP.show()