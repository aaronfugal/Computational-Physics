import numpy as np
import matplotlib.pyplot as plt

def LinearInterp(x, x1, y1, x2, y2):
    return y1 + ((y2-y1)/(x2-x1))*(x-x1)

def Find_CornHole_Derivs(s, t):
    # This function accepts a state (s) and time (t)
    # The state input breaks down into:
    # s[0] is position (r)
    # s[1] is velocity (v)
    # The function returns the state derivatives
    #   return[0] is dr/dt (v)
    #   return[1] is dv/dt (a)
    v_now = s[1]
    a_now = -k * np.linalg.norm(v_now)*v_now + a0
    return np.array([v_now, a_now])

def ode_step(s, dt, Find_Derivs, method):
    # s is the state coming in (in this case: [r, v])
    if(method == 'cromer'):
        # Cromer breaks into the routine (not all of the past is on the left-hand of the equation)
        # So there isn't a clean step here:
        r, v = s[0], s[1]
        Derivs = Find_Derivs(s, t)  # derivs[0] is v, derivs[1] is a
        v_new = v + Derivs[1]*dt
        r_new = r + v_new*dt
        s_new = np.array([r_new, v_new])
    elif(method == 'rk2'):
        k1 = dt * Find_Derivs(s, t)
        k2 = dt * Find_Derivs(s + k1, t + dt)
        s_new = s + 0.5*(k1 + k2)
    elif(method == 'rk4'):
        k1 = dt * Find_Derivs(s, t)
        k2 = dt * Find_Derivs(s + 0.5*k1, t + 0.5*dt)
        k3 = dt * Find_Derivs(s + 0.5*k2, t + 0.5*dt)
        k4 = dt * Find_Derivs(s + k3, t + dt)
        s_new = s + (k1 + 2*k2 + 2*k3 + k4)/6.0
    else:  # default to euler
        s_new = s + dt * Find_Derivs(s, t)
        
    return s_new

# Cornhole target specifics
cornhole_distance = 9.144  # (m) Pitch distance (30 ft)
cornhole_radius = 0.1524/2  # (m) 6/2 inch hole
cornhole_height = 0.2667  # (m) 10.5 inches - The angle of the board is 9.6 degrees
cornhole_target = np.array([cornhole_distance, cornhole_height])
target_ll = np.array([cornhole_distance - cornhole_radius*np.cos(9.6*np.pi/180.0), cornhole_height - cornhole_radius*np.sin(9.6*np.pi/180.0)])
target_ur = np.array([cornhole_distance + cornhole_radius*np.cos(9.6*np.pi/180.0), cornhole_height + cornhole_radius*np.sin(9.6*np.pi/180.0)])
mc = (target_ur[1] - target_ll[1])/(target_ur[0] - target_ll[0])

# Cornhole bag specifics
r_object = 0.1524  # m - Cornhole bean bags are 6in square and they are 
m_object = 0.449  # (kg) 15.85 oz bean bag

DiffEQMethod = 'rk4'  # Choices are 'rk4', 'rk2', 'cromer', 'euler (default)'

v_initial = 7.45  # m/s
theta = 18.0 # degrees
g = 3.8  # m/s^2
dt = 0.1  # seconds

print(f'\nInitial velocity {v_initial:0.1f} m/s at {theta:0.1f} degrees.')

a0 = np.array([0.0, -g])  # m/s^2
v0 = v_initial * np.array([np.cos(theta * np.pi/180.0), np.sin(theta * np.pi/180.0)])  # m/s
r0 = np.array([0.0, 0.5])  # m

total_time = 2 * -v0[1] / a0[1]
print(f'The time without air resistance for this trajectory is approximately {total_time:0.2f} seconds.')

est_time = int(total_time + 1)

t_start = 0.0  # seconds
t_end = est_time  # seconds

N = int((t_end-t_start)/dt + 1)
print(f'This will simulate from {t_start:0.1f} to {t_end:0.1f} seconds with a time spacing of {dt:0.4f} seconds.')
print(f'This will require {N} bins (including 0).')

C = 0.5  # Drag coefficient
rho_air = 0.020  # kg/m^3
k = (C * rho_air * np.pi * r_object**2) / (2 * m_object)
# k = 0.0 # No air resistance

t, v, r = np.linspace(t_start, t_end, N), np.zeros((N, 2)), np.zeros((N, 2))

v[0] = v0
r[0] = r0

i_uground = 0
i_won = 0

for i in range(1, N):
    s_now = np.array([r[i-1], v[i-1]])
    s_new = ode_step(s_now, dt, Find_CornHole_Derivs, DiffEQMethod)
    r[i], v[i] = s_new[0], s_new[1]
            
    # This checks to see if two bins straddle the line connecting the edges of the target hole
    if(r[i,0] > target_ll[0] and r[i,1] < target_ur[1]):  # Only bother if you are near the hole
        ms = (r[i,1] - r[i-1,1])/(r[i,0] - r[i-1,0])
        xi = (mc*target_ll[0] - ms*r[i-1,0] - target_ll[1] + r[i-1,1])/(mc-ms)
        yi = mc*(xi - target_ll[0]) + target_ll[1]
        # But count those only in the confines of the hole 
        # (don't count those that straddle the line away from the edges)
        if xi > target_ll[0] and xi < target_ur[0]:  # See if the projected crossing is near the hole
            if r[i-1, 0] < xi and r[i, 0] > xi:  # Make sure that the bins straddle the projected crossing
                i_won = i
                i_won_effective = LinearInterp(xi, r[i_won-1, 0], i_won-1, r[i_won, 0], i_won)
    
    # Check to see if we hit the ground, if so, quit and print where it landed
    if(r[i,1] < 0.0):
        i_uground = i
        print('We have hit the ground at around {r[i,0]:0.2f} m (between bins {i_uground-1} and {i_uground}).')
        break
    
if i_won != 0:
    print(f'I won at {t[i_won]:0.4f} seconds after releasing the throw (bin {i_won})!') 
    print(f'The effective bin of the win is {i_won_effective:0.4f} leading to a time of {i_won_effective*dt:0.4f} seconds.')
    print(f'The velocity in it\'s component form when it hits bin {i_won} is {v[14]} and the sum of the components has a')
    print(f'magnitude of 7.646 and an angle of 23.22 degrees below the x-axis')

    v_real = np.column_stack((v0[0] + a0[0]*t, v0[1] + a0[1]*t))
r_real = np.column_stack((r0[0] + v0[0]*t + 0.5*a0[0]*t**2, r0[1] + v0[1]*t + 0.5*a0[1]*t**2))

plt.figure('Position')
plt.title('Position')
plt.plot(r[0:i_uground+1,0], r[0:i_uground+1,1], 'b-', cornhole_distance, cornhole_height, 'r+')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()

from scipy.integrate import odeint

def Find_CornHole_Derivs_for_Solver(s, t):
    # This function accepts a state (s) and time (t)
    # The state input breaks down into:
    # s[0] is position-x (r_x)
    # s[1] is position-y (r_y)
    # s[2] is velocity-x (v_x)
    # s[3] is velocity-y (v_y)
    # The function returns the state derivatives
    #   return[0] is dr/dt-x (v_x)
    #   return[1] is dr/dt-y (v_y)
    #   return[2] is dv/dt-x (a_x)
    #   return[3] is dv/dt-y (a_y)
    vx_now, vy_now = s[2], s[3]
    v_now = np.array([vx_now, vy_now])
    ax_now = -k * np.linalg.norm(v_now)*vx_now + a0[0]
    ay_now = -k * np.linalg.norm(v_now)*vy_now + a0[1]
    return np.array([vx_now, vy_now, ax_now, ay_now])


s0 = np.array([r0[0], r0[1], v0[0], v0[1]])  # Initial state
s_solver = odeint(Find_CornHole_Derivs_for_Solver, s0, t)

r_solver = np.vstack([s_solver[:, 0], s_solver[:, 1]]).T
v_solver = np.vstack([s_solver[:, 2], s_solver[:, 3]]).T

plt.figure('Position')
plt.title('Position (odeint)')
plt.plot(r_solver[0:i_uground+1,0], r_solver[0:i_uground+1,1], 'b-', r_real[0:i_uground+1,0], r_real[0:i_uground+1,1], 'r--', cornhole_distance, cornhole_height, 'r+')
plt.xlabel('x (m)')
plt.ylabel('y (m)')
plt.show()