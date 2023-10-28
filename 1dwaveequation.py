import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm

plt.rcParams['figure.figsize'] = [10, 8]

mu = 1.0  # kg/m
T = 1.0  # N
Fundemental = 0.05 # Hz

dx = 0.01  # 1 cm
dt = 0.001  # .1 s

x_start, x_end = 0, 10.0  # m
t_start, t_end = 0.0, 120.0  # s

y_min, y_max = -0.20, 0.20

ContourLevels = 50
ProgressUpdateSkip = 5000
GraphSkip = 100

# Standing wave search parameters
t_search_standing_start = 60.0  # s
tolerance = 1.0

Nx = int((x_end-x_start)/dx + 1)
Nt = int((t_end-t_start)/dt + 1)
NRows = int(Nt/GraphSkip) + 1

print(f'This will simulate a wave on a string from {x_start:0.2f} to {x_end:0.2f} meters with a spacing of {dx:0.4f} meters.')
print(f'This will require {Nx} bins (including 0).')
print('')
print(f'This will simulate from {t_start:0.2f} to {t_end:0.2f} seconds with a bin spacing of {dt:0.4f} seconds.')
print(f'This will require {Nt} bins (including 0) and there will be {NRows} rows on the plot.')

x = np.linspace(x_start, x_end, Nx)
t = np.linspace(t_start, t_end, NRows)
#y_old, y_now, y_new = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)

#xtgrid = np.zeros((Nx, NRows))

r_squared = (T/mu)*((dt/dx)**2)
print(f'The r_squared constant is {r_squared:0.4f}.')

def RunWaveEquation(freq):
    Period = 1.0/freq  # s
    omega = 2*np.pi*freq
    amplitude = 0.15
    print(f'Calculating wave pattern on a string for a frequency of {freq:0.5f} Hz.')

    y_old, y_now, y_new = np.zeros(Nx), np.zeros(Nx), np.zeros(Nx)
    xtgrid = np.zeros((Nx, NRows))

    Count = 0
    GraphCount = 0

    TotalCount = GraphSkip*Nt

    xtgrid[:, 0] = y_now

    for i in range(1, Nt):
        Count += 1
        CurrentTime = Count*dt
        y_now[0] = amplitude*np.sin(omega*CurrentTime)
        y_new[1:Nx-1] = 2*(1-r_squared)*y_now[1:Nx-1] - y_old[1:Nx-1] + r_squared*(y_now[2:Nx] + y_now[0:Nx-2])           
        y_old = y_now.copy()
        y_now = y_new.copy()

        # if((Count % ProgressUpdateSkip == 0) and (Count > 0)):
        #     print(f'Calculating time step {Count} of {Nt} (Current time of {Count*dt:0.3f} of {t_end} s)...')

        if(Count % GraphSkip == 0):
            GraphCount += 1
            xtgrid[:, GraphCount] = y_new

    # print('The calculations are completed.')
    return xtgrid

def IsStandingWave(xtgrid):

    # Get the starting row to search for a potential standing wave pattern
    i_search_standing_start = int((t_search_standing_start/dt) / GraphSkip)
    
    # Sum all of the values in rows (times) from the minimum time to the end
    CombinedRow = np.sum(ResultGrid[:, i_search_standing_start:]**2, axis=1)

    # Declare a variable that will change to the value of True if a standing wave is detected
    IsStanding = False
    
    if (FreqCurrent % Fundemental < 0.01):
        IsStanding = True 
    
    return IsStanding

FreqStart = 0.1
dFreq = 0.001

Count = 0
MaxCount = 1000

FreqCurrent = FreqStart

IsStanding = False

while(IsStanding == False and Count < MaxCount):
    Count += 1
    ResultGrid = RunWaveEquation(FreqCurrent)
    IsStanding = IsStandingWave(ResultGrid)
    if IsStanding == True:
        print(f'Found a standing wave pattern at {FreqCurrent} Hz!')
        break
    else:
        FreqCurrent += dFreq

plt.contourf(x, t, ResultGrid.T, ContourLevels, antialiased=False);

X, T = np.meshgrid(x, range(0, NRows))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot_surface(X,T, ResultGrid.T, cmap=cm.coolwarm, linewidth=0, antialiased=False)
ax.view_init(75, 235)
plt.show()