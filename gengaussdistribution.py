### outputs a gauss function wave ###

import numpy as np
import matplotlib.pyplot as plt

### setup the time domain ###

N      =  1024   # power of 2, for best results   
tmin   = -15.0
tmax   =  15.0
deltat = (tmax-tmin)/float(N)

### parameters for the function ###

sigm = 1.0   # width of the gaussian
sigm = 0.25  # width of the gaussian
mean = 0.0   # mean of the gaussian

### step through domain, print t and f(t) ###

for i in range(0,N):
    time = tmin + deltat*i
    gaut = (1/np.sqrt(2.0*np.pi))*np.exp(-0.5*time*time)
    print("%16.12f %16.12f" % (time,gaut))
    
