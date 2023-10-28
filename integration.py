import numpy as np
import random
import matplotlib.pyplot as plt

NTHROWS = 1000  # for Monte Carlo method

###########
def fn(xx):
###########
    
    return np.sin(xx)

###################
### main script ###
###################

xmin =  0.0
xmax =  2.0

n = int(input("enter number of bins: "))
if(n%2): n += 1  # fix user mistake, make n even
print("Using n =",n)
print("NTHROWS =",NTHROWS)

dx   = (xmax-xmin)/n

### integral by trapezoid rule ### 

Ftrap = 0

for i in range(0,n+1):
    x = xmin + i*dx
    w = 1.0
    if(i==0 or i==n): w = 0.5
    Ftrap += w * fn(x) * dx

### integral by midpoint rule ### 

Fmidp = 0

for i in range(0,n):
    x = xmin + (i+0.5)*dx
    Fmidp += fn(x) * dx

### integral by Simpson's rule ### 
###    (lazy man's version)    ###

Fsimp = (2.0*Fmidp + Ftrap)/3.0

### integral done by Monte Carlo ### 

xmc = np.zeros(NTHROWS)
ymc = np.zeros(NTHROWS)

nunder = 0

for ithrow in range(0,NTHROWS):
    xmc[ithrow] = random.uniform(xmin,xmax)
    ymc[ithrow] = random.uniform(0.0,1.0)
    if(ymc[ithrow] < fn(xmc[ithrow])):
        nunder += 1

Fmont = (nunder/NTHROWS)*(xmax-xmin)
        
### plot the function being integrated (for debugging) ###

if(False):
    xplot = np.arange(xmin,xmax+dx,dx)
    yplot = fn(xplot)
    plt.plot(xplot,yplot)
    plt.scatter(xmc,ymc,c='k',s=0.2)
    plt.grid()
    plt.xlabel("$x$")
    plt.ylabel("$f(x)$")
    plt.show()

### integral done analytically ### 

Fanal = np.cos(xmin) - np.cos(xmax)

### output results to terminal ### 

print("Integral using trapezoid rule = ",Ftrap)
print("Integral using midpoint rule  = ",Fmidp)
print("Integral using Simpson's rule = ",Fsimp)
print("Integral using Monte Carlo    = ",Fmont)
print("Integral solved analytically  = ",Fanal)


#