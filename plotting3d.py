import numpy as np
from mpl_toolkits.mplot3d import Axes3D  
import matplotlib.pyplot as plt

### the 2D function to be plotted ###
def fun(x, y):
    return x**2 + y

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
x = y = np.arange(-3.0, 3.0, 0.05)

# meshgrid produces 2D arrays X and Y, whose elements 
# are given the components of vectors x and y respectively
X, Y = np.meshgrid(x, y)

Z = fun(X,Y)

ax.plot_surface(X, Y, Z)

ax.set_xlabel('X Label')
ax.set_ylabel('Y Label')
ax.set_zlabel('Z Label')

plt.show()