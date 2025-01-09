import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

#constants/cond
g = 9.81  # acceleration due to gravity (m/s^2)


x0 = 0      # initial x position (m)
y0 = 0      # initial y position (m)
vx0 = 10    # initial x velocity (m/s)
vy0 = 10    # initial y velocity (m/s)
initial_state = [x0, y0, vx0, vy0]

# time para
t0 = 0
tf = 2 * vy0 / g  # time of flight
n_steps = 100
dt = (tf - t0) / n_steps


t_values = np.linspace(t0, tf, n_steps + 1)
x_values = np.zeros(n_steps + 1)
y_values = np.zeros(n_steps + 1)
vx_values = np.zeros(n_steps + 1)
vy_values = np.zeros(n_steps + 1)


x_values[0], y_values[0] = x0, y0
vx_values[0], vy_values[0] = vx0, vy0

# euler method 
for i in range(n_steps):
    x_values[i + 1] = x_values[i] + dt * vx_values[i]
    y_values[i + 1] = y_values[i] + dt * vy_values[i]
    vx_values[i + 1] = vx_values[i]
    vy_values[i + 1] = vy_values[i] #- dt * g


fig, ax = plt.subplots()
ax.set_xlim(0, max(x_values) + 5)
ax.set_ylim(0, max(y_values) + 5)
point, = ax.plot(x_values[0], y_values[0], 'bo')

def init():
    point.set_data(x_values[0], y_values[0])
    return point,

def update(frame):
    point.set_data(x_values[frame], y_values[frame])
    return point,

ani = FuncAnimation(fig, update, frames=len(x_values), init_func=init, blit=True, interval=25)

plt.show()
