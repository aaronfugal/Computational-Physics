import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

G = 6.67e-11  # Gravitational Constant [Nm^2/kg^2]

# Earth Information
R_earth = 6.37e6  # Radius of the Earth [m]
M_earth = 5.972e24  # Mass of the Earth [kg]


# Moon Information

R_moon = 1.74e6  # Radius of the Moon [m]
M_moon = 7.36e22  # Mass of the moon [kg]

# Script Information
ContourLevels = 100  # This can be much higher if needed (up to 1000)
ProgressUpdateSkip = 100

x_start, y_start = 0.0, -400e5  # m
x_stop, y_stop = 400e6, 400e6  # m

h = 1000000

Nx, Ny = int((x_stop-x_start)/h + 1), int((y_stop-y_start)/h + 1)
print(f'This will simulate from x = {x_start} to {x_stop} distance units with a spacing of {h:0.4f} units.')
print(f'This will simulate from y = {y_start} to {y_stop} distance units with a spacing of {h:0.4f} units.')
print(f'This will require {Nx} by {Ny} bins (including 0).')

def EarthDist(x, y):
    return(np.sqrt(x**2 + y**2))

def MoonDist(x, y):
    return(np.sqrt((x-8e7)**2 + y**2))

x, y = np.linspace(x_start, x_stop, Nx), np.linspace(y_start, y_stop, Ny)
xygrid = np.zeros((Nx, Ny))

RowCount = 0

for i in range(Nx):
    RowCount += 1
    if(RowCount % ProgressUpdateSkip == 0):
        print(f'Currently on iteration {RowCount} of {Nx}...')
    for j in range(Ny):
        if(EarthDist(x[i], y[j]) >= R_earth) and (MoonDist(x[i], y[j]) >= R_moon):
            xygrid[i, j] = -G*M_earth/EarthDist(x[i], y[j]) + -G*M_moon/MoonDist(x[i], y[j])
        else:
            xygrid[i, j] = 0.0

print(f'\nFinished.')
plt.contourf(x, y, xygrid.T, ContourLevels, cmap=cm.jet);
plt.show()
plt.contour(x, y, xygrid.T, ContourLevels, cmap=cm.jet);
plt.show()

