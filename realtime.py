
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


#const/cond
#increasing case: throguh the bottom, decreasing case through the top, show for each pt


G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_jupiter = 1898e24  # mass of the jupiter, kg
M_satellite = 1000  # mass of the satellite, kg (arbitrarily chosen for visualization)
R_jupiter = 71492e3  # radius of the jupiter, m


r0 = R_jupiter + 1.5e10  # initial satellite distance from jupiter's center, m 
ve0 = 13.1e3 * 0  # m/s, jupiter's initial speed
vs0 = 10e3 + 15e3  # earth speed +  initial speed for a circular orbit, m/s (e.g. Pinoeer at 56 km/s, fastest)
stheta = -35 - .8 # angle of satellite compared to Jupiter (degrees)


initial_state = [0, 0.5e10, 0, 0, -r0, r0*1.1, vs0*np.cos(np.radians(stheta)), vs0*np.sin(np.radians(stheta))]  # [x_jupiter, y_jupiter, vx_jupiter, vy_jupiter, x_sat, y_sat, vx_sat, vy_sat]

# time para
t0 = 0
tf = (4*r0)/vs0
n_steps = 50000
dt = (tf - t0) / n_steps

t_values = np.linspace(t0, tf, n_steps + 1)
x_jupiter_values = np.zeros(n_steps + 1)
y_jupiter_values = np.zeros(n_steps + 1)
vx_jupiter_values = np.zeros(n_steps + 1)
vy_jupiter_values = np.zeros(n_steps + 1)
x_sat_values = np.zeros(n_steps + 1)
y_sat_values = np.zeros(n_steps + 1)
vx_sat_values = np.zeros(n_steps + 1)
vy_sat_values = np.zeros(n_steps + 1)

x_jupiter_values[0], y_jupiter_values[0] = initial_state[0], initial_state[1]
vx_jupiter_values[0], vy_jupiter_values[0] = initial_state[2], initial_state[3]
x_sat_values[0], y_sat_values[0] = initial_state[4], initial_state[5]
vx_sat_values[0], vy_sat_values[0] = initial_state[6], initial_state[7]

# velocity verlet (newton grav)
# ax = g * m * (delta x) / r^3 
# ay = g * m * (delta y) / r^3 
# x = x + vx * timestep + 0.5 * ax * time step^2 
# y = y + yx * timestep + 0.5 * ay * time step^2 
# delta v = 0.5 * (a[t] + a[t+1]) * time step 
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
    input_file.write(f"{t_values[i]} {np.sqrt(vx_sat_values[i]**2+ vy_sat_values[i]**2)} {vx_sat_values[i]} {vy_sat_values[i]} {x_sat_values[i]} {y_sat_values[i]}\n")

input_file.close()


fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))
ax1.set_aspect('equal')
scale_lim = 10
ax1.set_xlim(-scale_lim*2.5e9, scale_lim*2.5e9)
ax1.set_ylim(1e9, scale_lim*3e9)
jupiter, = ax1.plot([], [], 'o', color='b', markersize=10) 
satellite, = ax1.plot([], [], 'o', color='red')
trail, = ax1.plot([], [], '-', color='Black', lw=1, alpha=0.5)

ax2.set_xlim(-scale_lim*r0/10, scale_lim*r0/10)
ax2.set_ylim(20000, np.max(np.sqrt(vx_sat_values**2 + vy_sat_values**2)) + 5000)

speed_line, = ax2.plot([], [], 'r-', label='Speed (m/s)')
ax2.axhline(y=25000, color='blue', linestyle='--', label='Base Speed') 

ax2.set_xlabel('Position x (m)')
ax2.set_ylabel('Speed (m/s)')
ax2.legend()





def init():
    jupiter.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    speed_line.set_data([], [])
    return jupiter, satellite, trail, speed_line

def update(frame):
    frame = (frame * 100) % n_steps
    jupiter.set_data(x_jupiter_values[frame], y_jupiter_values[frame])
    satellite.set_data(x_sat_values[frame], y_sat_values[frame])
    trail.set_data(x_sat_values[:frame], y_sat_values[:frame])
   
    speed = np.sqrt(vx_sat_values[frame]**2 + vy_sat_values[frame]**2)
    speed_line.set_data(x_sat_values[:frame], np.sqrt(vx_sat_values[:frame]**2 + vy_sat_values[:frame]**2))
   
    return jupiter, satellite, trail, speed_line

ani = FuncAnimation(fig, update, frames=len(t_values), init_func=init, blit=True, interval=0)

print(f"Starting speed:{np.sqrt(vx_sat_values[0]**2 + vy_sat_values[0]**2)}, Ending speed:{np.sqrt(vx_sat_values[n_steps]**2 + vy_sat_values[n_steps]**2)}")

plt.show()



