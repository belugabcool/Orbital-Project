import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# const/cond
k = 10    # spring constant (N/m)
m = 1     # mass of the object (kg)

length_equilibrium = 1


x0 = length_equilibrium + 0.2  # initial position (m)
v0 = 0                         # initial velocity (m/s)
initial_state = [x0, v0]

# time par 
t0 = 0
tf = 10  # 10 seconds
n_steps = 300
dt = (tf - t0) / n_steps


t = np.linspace(t0, tf, n_steps + 1)
x = np.zeros(n_steps + 1)
v = np.zeros(n_steps + 1)


x[0], v[0] = x0, v0

# velocity verlet
for i in range(n_steps):
    a = -(k/m) * (x[i] - length_equilibrium)
    x[i + 1] = x[i] + v[i] * dt + 0.5 * a * dt**2
    a_new = -(k/m) * (x[i + 1] - length_equilibrium)
    v[i + 1] = v[i] + 0.5 * (a + a_new) * dt

fig, ax = plt.subplots()
ax.set_xlim(length_equilibrium - 0.5, length_equilibrium + 0.5)
ax.set_ylim(-1, 1)
line, = ax.plot([length_equilibrium, x[0]], [0, 0], color='gray')  # Spring line
mass, = ax.plot(x[0], 0, 'bo', markersize=10)  # Mass

def init():
    line.set_data([length_equilibrium, x[0]], [0, 0])
    mass.set_data(x[0], 0)
    return line, mass

def update(frame):
    line.set_data([length_equilibrium, x[frame]], [0, 0])
    mass.set_data(x[frame], 0)
    return line, mass

ani = FuncAnimation(fig, update, frames=len(x), init_func=init, blit=True, interval=50)

plt.show()
