
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


# Constants


G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_jupiter = 1898e24  # mass of the jupiter, kg
M_satellite = 1000  # mass of the satellite, kg (arbitrarily chosen for visualization)
R_jupiter = 71492e3  # radius of the jupiter, m

# Initial conditions
r0 = R_jupiter + 1.5e10  # initial satellite distance from jupiter's center, m 
ve0 = 13.1e3  # m/s, jupiter's initial speed
vs0 = 40e3  # initial speed for a circular orbit, m/s
stheta = 12 # angle of satellite compared to Jupiter (degrees)


initial_state = [0, 0, 0, ve0, -r0, 0, vs0*np.cos(np.radians(stheta)), vs0*np.sin(np.radians(stheta))]  # [x_jupiter, y_jupiter, vx_jupiter, vy_jupiter, x_sat, y_sat, vx_sat, vy_sat]

# Time parameters
t0 = 0
tf = (2*r0)/vs0
n_steps = 3000
dt = (tf - t0) / n_steps

# Initialize arrays to store the solution
t_values = np.linspace(t0, tf, n_steps + 1)
x_jupiter_values = np.zeros(n_steps + 1)
y_jupiter_values = np.zeros(n_steps + 1)
vx_jupiter_values = np.zeros(n_steps + 1)
vy_jupiter_values = np.zeros(n_steps + 1)
x_sat_values = np.zeros(n_steps + 1)
y_sat_values = np.zeros(n_steps + 1)
vx_sat_values = np.zeros(n_steps + 1)
vy_sat_values = np.zeros(n_steps + 1)

# Set initial conditions
x_jupiter_values[0], y_jupiter_values[0] = initial_state[0], initial_state[1]
vx_jupiter_values[0], vy_jupiter_values[0] = initial_state[2], initial_state[3]
x_sat_values[0], y_sat_values[0] = initial_state[4], initial_state[5]
vx_sat_values[0], vy_sat_values[0] = initial_state[6], initial_state[7]


for i in range(n_steps):
    r = np.sqrt((x_sat_values[i] - x_jupiter_values[i])**2 + (y_sat_values[i] - y_jupiter_values[i])**2)
    
    dist = np.sqrt((x_sat_values[i]-x_jupiter_values[i])**2 + (y_sat_values[i]-y_jupiter_values[i])**2)
    if dist < R_jupiter:
        break 
    
    ax_jupiter = G * M_satellite * (x_sat_values[i] - x_jupiter_values[i]) / r**3
    ay_jupiter = G * M_satellite * (y_sat_values[i] - y_jupiter_values[i]) / r**3
    ax_sat = -G * M_jupiter * (x_sat_values[i] - x_jupiter_values[i]) / r**3
    ay_sat = -G * M_jupiter * (y_sat_values[i] - y_jupiter_values[i]) / r**3
    
    x_jupiter_values[i + 1] = x_jupiter_values[i] + vx_jupiter_values[i] * dt + 0.5 * ax_jupiter * dt**2
    y_jupiter_values[i + 1] = y_jupiter_values[i] + vy_jupiter_values[i] * dt + 0.5 * ay_jupiter * dt**2
    x_sat_values[i + 1] = x_sat_values[i] + vx_sat_values[i] * dt + 0.5 * ax_sat * dt**2
    y_sat_values[i + 1] = y_sat_values[i] + vy_sat_values[i] * dt + 0.5 * ay_sat * dt**2
    
    r_new = np.sqrt((x_sat_values[i + 1] - x_jupiter_values[i + 1])**2 + (y_sat_values[i + 1] - y_jupiter_values[i + 1])**2)
    ax_jupiter_new = G * M_satellite * (x_sat_values[i + 1] - x_jupiter_values[i + 1]) / r_new**3
    ay_jupiter_new = G * M_satellite * (y_sat_values[i + 1] - y_jupiter_values[i + 1]) / r_new**3
    ax_sat_new = -G * M_jupiter * (x_sat_values[i + 1] - x_jupiter_values[i + 1]) / r_new**3
    ay_sat_new = -G * M_jupiter * (y_sat_values[i + 1] - y_jupiter_values[i + 1]) / r_new**3
    

    vx_jupiter_values[i + 1] = vx_jupiter_values[i] + 0.5 * (ax_jupiter + ax_jupiter_new) * dt
    vy_jupiter_values[i + 1] = vy_jupiter_values[i] + 0.5 * (ay_jupiter + ay_jupiter_new) * dt
    vx_sat_values[i + 1] = vx_sat_values[i] + 0.5 * (ax_sat + ax_sat_new) * dt
    vy_sat_values[i + 1] = vy_sat_values[i] + 0.5 * (ay_sat + ay_sat_new) * dt
    


input_file = open("dat.txt", "w")
for i in range(n_steps):
    input_file.write(f"{t_values[i]} {vx_sat_values[i]} {vy_sat_values[i]} {x_sat_values[i]} {y_sat_values[i]}\n")

input_file.close()

fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
scale_lim = 1.5
ax.set_xlim(-scale_lim*r0, scale_lim*r0)
ax.set_ylim(-scale_lim*r0, scale_lim*r0)
jupiter, = ax.plot([], [], 'o', color='b', markersize=10)  
satellite, = ax.plot([], [], 'o', color='red')
trail, = ax.plot([], [], '-', color='green', lw=1, alpha=0.5)

def init():
    jupiter.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    return jupiter, satellite, trail

def update(frame):
    jupiter.set_data(x_jupiter_values[frame], y_jupiter_values[frame])
    satellite.set_data(x_sat_values[frame], y_sat_values[frame])
    trail.set_data(x_sat_values[:frame], y_sat_values[:frame])
    return jupiter, satellite, trail

ani = FuncAnimation(fig, update, frames=len(t_values), init_func=init, blit=True, interval=1)

plt.show()
