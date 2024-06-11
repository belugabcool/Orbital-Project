import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Constants
k = 50    # spring constant (N/m)
m = 1     # mass of the object (kg)

# Equilibrium length of the spring (m)
length_equilibrium = 1

# Initial conditions
# Assume the object is initially stretched 0.2 meters to the right of equilibrium and released with no initial velocity
x0 = length_equilibrium + 1  # initial position (m)
v0 = 0                       # initial velocity (m/s)
initial_state = [x0, v0]

# Time parameters
t0 = 0
tf = 10  # 10 seconds
n_steps = 300
dt = (tf - t0) / n_steps

# Initialize arrays to store the solution
t_values = np.linspace(t0, tf, n_steps + 1)
x_values = np.zeros(n_steps + 1)
v_values = np.zeros(n_steps + 1)

# Set initial conditions
x_values[0], v_values[0] = x0, v0

# Euler method to solve the equations of motion
for i in range(n_steps):
    a = -(k/m) * (x_values[i] - length_equilibrium)
    v_values[i + 1] = v_values[i] + dt * a
    x_values[i + 1] = x_values[i] + dt * v_values[i]
    

# Create the animation
fig, ax = plt.subplots()
ax.set_xlim(length_equilibrium - 1, length_equilibrium + 1)
ax.set_ylim(-1, 1)
line, = ax.plot([length_equilibrium, x_values[0]], [0, 0], color='gray')  # Spring line
mass, = ax.plot(x_values[0], 0, 'bo', markersize=10)  # Mass

def init():
    line.set_data([length_equilibrium, x_values[0]], [0, 0])
    mass.set_data(x_values[0], 0)
    return line, mass

def update(frame):
    line.set_data([length_equilibrium, x_values[frame]], [0, 0])
    mass.set_data(x_values[frame], 0)
    return line, mass

ani = FuncAnimation(fig, update, frames=len(x_values), init_func=init, blit=True, interval=50)

plt.show()
