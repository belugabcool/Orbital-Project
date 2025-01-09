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


t = np.linspace(t0, tf, n_steps + 1)
x = np.zeros(n_steps + 1)
y = np.zeros(n_steps + 1)
vx = np.zeros(n_steps + 1)
vy = np.zeros(n_steps + 1)


x[0], y[0] = x0, y0
vx[0], vy[0] = vx0, vy0

# euler method 
for i in range(n_steps):
    x[i + 1] = x[i] + dt * vx[i]
    y[i + 1] = y[i] + dt * vy[i]
    vx[i + 1] = vx[i]
    vy[i + 1] = vy[i] #- dt * g


fig, ax = plt.subplots()
ax.set_xlim(0, max(x) + 5)
ax.set_ylim(0, max(y) + 5)
point, = ax.plot(x[0], y[0], 'bo')

def init():
    point.set_data(x[0], y[0])
    return point,

def update(frame):
    point.set_data(x[frame], y[frame])
    return point,

ani = FuncAnimation(fig, update, frames=len(x), init_func=init, blit=True, interval=25)

plt.show()
