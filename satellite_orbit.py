import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#constants/conditions
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M = 5.972e24  # mass of the earth, kg
R = 6.371e6  # radius of the earth, m

r0 = R + 2000e3 # initial satellite distance from Earth's center, m (400 km altitude)
# v0 = np.sqrt(G * M / r0) * 1.0  # initial speed for a near-circular orbit, m/s
v0 =  8000  # initial speed for a near-circular orbit, m/s
initial_state = [r0, 0, 0, v0]  # [x0, y0, vx0, vy0]
frameScale = 7

# time vals
t0 = 0
orbit_period = 10 * 2 * np.pi * np.sqrt(r0**3 / (G * M)) 
tf = orbit_period * 10
n_steps = 3000
dt = (tf - t0) / n_steps

t_v = np.linspace(t0, tf, n_steps + 1)
x_v = np.zeros(n_steps + 1)
y_v = np.zeros(n_steps + 1)
vx_v = np.zeros(n_steps + 1)
vy_v = np.zeros(n_steps + 1)

x_v[0], y_v[0] = initial_state[0], initial_state[1]
vx_v[0], vy_v[0] = initial_state[2], initial_state[3]

# velocity verlet
for i in range(n_steps):
    r = np.sqrt(x_v[i]**2 + y_v[i]**2)
    ax = -G * M * x_v[i] / r**3
    ay = -G * M * y_v[i] / r**3
    
    x_v[i + 1] = x_v[i] + vx_v[i] * dt + 0.5 * ax * dt**2
    y_v[i + 1] = y_v[i] +vy_v[i]* dt + 0.5 * ay * dt**2
    
    r_new = np.sqrt(x_v[i + 1]**2 + y_v[i + 1]**2)
    ax_new = -G * M * x_v[i + 1] / r_new**3
    ay_new = -G * M * y_v[i + 1] / r_new**3
    
    vx_v[i + 1] = vx_v[i] + 0.5 * (ax + ax_new) * dt
    vy_v[i + 1] = vy_v[i] + 0.5 * (ay + ay_new) * dt

# animations
fig, ax = plt.subplots(figsize=(8, 8))
ax.set_aspect('equal')


ax.set_xlim(-frameScale*r0, frameScale*r0)
ax.set_ylim(-frameScale*r0, frameScale*r0)
earth = plt.Circle((0, 0), R, color='b')
satellite, = ax.plot([], [], 'o', color='red')
trail, = ax.plot([], [], '-', color='yellow', lw=1, alpha=0.5)
ax.add_artist(earth)

def init():
    satellite.set_data([], [])
    trail.set_data([], [])
    return satellite, trail

def update(frame):
    satellite.set_data(x_v[frame], y_v[frame])
    trail.set_data(x_v[:frame], y_v[:frame])
    return satellite, trail

ani = FuncAnimation(fig, update, frames=len(t_v), init_func=init, blit=True, interval=10)

plt.show()
