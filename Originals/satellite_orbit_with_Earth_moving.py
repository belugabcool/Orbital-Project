import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_earth = 5.972e24  # mass of the Earth, kg
M_satellite = 500  # mass of the satellite, kg (arbitrarily chosen for visualization)
R_earth = 6.371e6  # radius of the Earth, m

# Initial conditions
r0 = R_earth + 400e3  # initial satellite distance from Earth's center, m (400 km altitude)
ve0 = 1000  # m/s, Earth's initial speed
vs0 = ve0 + np.sqrt(G * M_earth / r0)  # initial speed for a circular orbit, m/s
initial_state = [0, 0, 0, ve0, r0, 0, 0, vs0]  # [x_earth, y_earth, vx_earth, vy_earth, x_sat, y_sat, vx_sat, vy_sat]

# Time parameters
t0 = 0
orbit_period = 10 * 2 * np.pi * np.sqrt(r0**3 / (G * M_earth))
tf = orbit_period
n_steps = 3000
dt = (tf - t0) / n_steps

# Initialize arrays to store the solution
t_values = np.linspace(t0, tf, n_steps + 1)
x_earth_values = np.zeros(n_steps + 1)
y_earth_values = np.zeros(n_steps + 1)
vx_earth_values = np.zeros(n_steps + 1)
vy_earth_values = np.zeros(n_steps + 1)
x_sat_values = np.zeros(n_steps + 1)
y_sat_values = np.zeros(n_steps + 1)
vx_sat_values = np.zeros(n_steps + 1)
vy_sat_values = np.zeros(n_steps + 1)

# Set initial conditions
x_earth_values[0], y_earth_values[0] = initial_state[0], initial_state[1]
vx_earth_values[0], vy_earth_values[0] = initial_state[2], initial_state[3]
x_sat_values[0], y_sat_values[0] = initial_state[4], initial_state[5]
vx_sat_values[0], vy_sat_values[0] = initial_state[6], initial_state[7]

# Velocity Verlet method to solve the equations of motion
for i in range(n_steps):
    r = np.sqrt((x_sat_values[i] - x_earth_values[i])**2 + (y_sat_values[i] - y_earth_values[i])**2)
    
    # Accelerations
    ax_earth = G * M_satellite * (x_sat_values[i] - x_earth_values[i]) / r**3
    ay_earth = G * M_satellite * (y_sat_values[i] - y_earth_values[i]) / r**3
    ax_sat = -G * M_earth * (x_sat_values[i] - x_earth_values[i]) / r**3
    ay_sat = -G * M_earth * (y_sat_values[i] - y_earth_values[i]) / r**3
    
    # Update positions
    x_earth_values[i + 1] = x_earth_values[i] + vx_earth_values[i] * dt + 0.5 * ax_earth * dt**2
    y_earth_values[i + 1] = y_earth_values[i] + vy_earth_values[i] * dt + 0.5 * ay_earth * dt**2
    x_sat_values[i + 1] = x_sat_values[i] + vx_sat_values[i] * dt + 0.5 * ax_sat * dt**2
    y_sat_values[i + 1] = y_sat_values[i] + vy_sat_values[i] * dt + 0.5 * ay_sat * dt**2
    
    # Calculate new accelerations
    r_new = np.sqrt((x_sat_values[i + 1] - x_earth_values[i + 1])**2 + (y_sat_values[i + 1] - y_earth_values[i + 1])**2)
    ax_earth_new = G * M_satellite * (x_sat_values[i + 1] - x_earth_values[i + 1]) / r_new**3
    ay_earth_new = G * M_satellite * (y_sat_values[i + 1] - y_earth_values[i + 1]) / r_new**3
    ax_sat_new = -G * M_earth * (x_sat_values[i + 1] - x_earth_values[i + 1]) / r_new**3
    ay_sat_new = -G * M_earth * (y_sat_values[i + 1] - y_earth_values[i + 1]) / r_new**3
    
    # Update velocities
    vx_earth_values[i + 1] = vx_earth_values[i] + 0.5 * (ax_earth + ax_earth_new) * dt
    vy_earth_values[i + 1] = vy_earth_values[i] + 0.5 * (ay_earth + ay_earth_new) * dt
    vx_sat_values[i + 1] = vx_sat_values[i] + 0.5 * (ax_sat + ax_sat_new) * dt
    vy_sat_values[i + 1] = vy_sat_values[i] + 0.5 * (ay_sat + ay_sat_new) * dt

# Create the animation
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')
scale_lim = 1.5
ax.set_xlim(-scale_lim*r0, scale_lim*r0)
ax.set_ylim(-scale_lim*r0, scale_lim*r0)
earth, = ax.plot([], [], 'o', color='b', markersize=10)  # Exaggerated size for visibility
satellite, = ax.plot([], [], 'o', color='red')
trail, = ax.plot([], [], '-', color='yellow', lw=1, alpha=0.5)

def init():
    earth.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    return earth, satellite, trail

def update(frame):
    earth.set_data(x_earth_values[frame], y_earth_values[frame])
    satellite.set_data(x_sat_values[frame], y_sat_values[frame])
    trail.set_data(x_sat_values[:frame], y_sat_values[:frame])
    return earth, satellite, trail

ani = FuncAnimation(fig, update, frames=len(t_values), init_func=init, blit=True, interval=50)

plt.show()
