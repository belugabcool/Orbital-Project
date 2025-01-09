import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# const/cond
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_earth = 5.972e24  # mass of the Earth, kg
M_satellite = 500  # mass of the satellite, kg (arbitrarily chosen for visualization)
R_earth = 6.371e6  # radius of the Earth, m


r0 = R_earth + 400e3  # initial satellite distance from Earth's center, m (400 km altitude)
ve0 = 1000  # m/s, Earth's initial speed
vs0 = ve0 + np.sqrt(G * M_earth / r0)  # initial speed for a circular orbit, m/s
initial_state = [0, 0, 0, ve0, r0, 0, 0, vs0]  # [x_earth, y_earth, vx_earth, vy_earth, x_sat, y_sat, vx_sat, vy_sat]

# time para
t0 = 0
orbit_period = 10 * 2 * np.pi * np.sqrt(r0**3 / (G * M_earth))
tf = orbit_period
n_steps = 3000
dt = (tf - t0) / n_steps


t = np.linspace(t0, tf, n_steps + 1)
x_earth = np.zeros(n_steps + 1)
y_earth = np.zeros(n_steps + 1)
vx_earth = np.zeros(n_steps + 1)
vy_earth = np.zeros(n_steps + 1)
x_sat = np.zeros(n_steps + 1)
y_sat = np.zeros(n_steps + 1)
vx_sat = np.zeros(n_steps + 1)
vy_sat = np.zeros(n_steps + 1)


x_earth[0], y_earth[0] = initial_state[0], initial_state[1]
vx_earth[0], vy_earth[0] = initial_state[2], initial_state[3]
x_sat[0], y_sat[0] = initial_state[4], initial_state[5]
vx_sat[0], vy_sat[0] = initial_state[6], initial_state[7]

# velocity verlet (newton grav)
# ax = g * m * (delta x) / r^3 
# ay = g * m * (delta y) / r^3 
# x = x + vx * timestep + 0.5 * ax * time step^2 
# y = y + yx * timestep + 0.5 * ay * time step^2 
# delta v = 0.5 * (a[t] + a[t+1]) * time step 
for i in range(n_steps):
    r = np.sqrt((x_sat[i] - x_earth[i])**2 + (y_sat[i] - y_earth[i])**2)
    
    ax_earth = G * M_satellite * (x_sat[i] - x_earth[i]) / r**3
    ay_earth = G * M_satellite * (y_sat[i] - y_earth[i]) / r**3
    ax_sat = -G * M_earth * (x_sat[i] - x_earth[i]) / r**3
    ay_sat = -G * M_earth * (y_sat[i] - y_earth[i]) / r**3
    
    x_earth[i + 1] = x_earth[i] + vx_earth[i] * dt + 0.5 * ax_earth * dt**2
    y_earth[i + 1] = y_earth[i] + vy_earth[i] * dt + 0.5 * ay_earth * dt**2
    x_sat[i + 1] = x_sat[i] + vx_sat[i] * dt + 0.5 * ax_sat * dt**2
    y_sat[i + 1] = y_sat[i] + vy_sat[i] * dt + 0.5 * ay_sat * dt**2
    
    r_new = np.sqrt((x_sat[i + 1] - x_earth[i + 1])**2 + (y_sat[i + 1] - y_earth[i + 1])**2)
    ax_earth_new = G * M_satellite * (x_sat[i + 1] - x_earth[i + 1]) / r_new**3
    ay_earth_new = G * M_satellite * (y_sat[i + 1] - y_earth[i + 1]) / r_new**3
    ax_sat_new = -G * M_earth * (x_sat[i + 1] - x_earth[i + 1]) / r_new**3
    ay_sat_new = -G * M_earth * (y_sat[i + 1] - y_earth[i + 1]) / r_new**3
    z
    vx_earth[i + 1] = vx_earth[i] + 0.5 * (ax_earth + ax_earth_new) * dt
    vy_earth[i + 1] = vy_earth[i] + 0.5 * (ay_earth + ay_earth_new) * dt
    vx_sat[i + 1] = vx_sat[i] + 0.5 * (ax_sat + ax_sat_new) * dt
    vy_sat[i + 1] = vy_sat[i] + 0.5 * (ay_sat + ay_sat_new) * dt

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
scale_lim = 1.5
ax.set_xlim(-scale_lim*r0, scale_lim*r0)
ax.set_ylim(-scale_lim*r0, scale_lim*r0)
earth, = ax.plot([], [], 'o', color='b', markersize=10)  
satellite, = ax.plot([], [], 'o', color='red')
trail, = ax.plot([], [], '-', color='yellow', lw=1, alpha=0.5)

def init():
    earth.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    return earth, satellite, trail

def update(frame):
    earth.set_data(x_earth[frame], y_earth[frame])
    satellite.set_data(x_sat[frame], y_sat[frame])
    trail.set_data(x_sat[:frame], y_sat[:frame])
    return earth, satellite, trail

ani = FuncAnimation(fig, update, frames=len(t), init_func=init, blit=True, interval=50)

plt.show()
