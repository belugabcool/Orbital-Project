
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from PIL import Image, ImageSequence

# Constants
#increasing case: throguh the bottom, decreasing case through the top, show for each pt
# D:\Users\lkhom\Desktop\PythonDS\Orbital-Project\boom2.gif

G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_jupiter = 1898e24  # mass of the jupiter, kg
M_satellite = 1000  # mass of the satellite, kg (arbitrarily chosen for visualization)
R_jupiter = 71492e3  # radius of the jupiter, m

# Initial conditions
r0 = R_jupiter + 1.5e10  # initial satellite distance from jupiter's center, m 
ve0 = 13.1e3 # m/s, jupiter's initial speed
vs0 = 40e3 # earth speed +  initial speed for a circular orbit, m/s (e.g. Pinoeer at 56 km/s, fastest)
stheta = 32.1 # angle of satellite compared to Jupiter (degrees)


initial_state = [0, 0, 0, ve0, -r0, -r0*0.25, vs0*np.cos(np.radians(stheta)), vs0*np.sin(np.radians(stheta))]  # [x_jupiter, y_jupiter, vx_jupiter, vy_jupiter, x_sat, y_sat, vx_sat, vy_sat]

# Time parameters
t0 = 0
tf = (5*r0)/vs0
n_steps = 50000
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

# Velocity Verlet method to solve the equations of motion
collision_check = False

for i in range(n_steps):
    r = np.sqrt((x_sat_values[i] - x_jupiter_values[i])**2 + (y_sat_values[i] - y_jupiter_values[i])**2)
    
    # check for collision between satellite and jupiter
    if r < R_jupiter:
        collision_check = True

    if(collision_check):
        vx_jupiter_values[i + 1] = 0
        vy_jupiter_values[i + 1] = 0
        vx_sat_values[i + 1] = 0
        vy_sat_values[i + 1] = 0

        x_jupiter_values[i + 1] = x_jupiter_values[i]
        y_jupiter_values[i + 1] = y_jupiter_values[i] 
        x_sat_values[i + 1] = x_sat_values[i]
        y_sat_values[i + 1] = y_sat_values[i]

        continue 
    


    # Accelerations
    ax_jupiter = G * M_satellite * (x_sat_values[i] - x_jupiter_values[i]) / r**3
    ay_jupiter = G * M_satellite * (y_sat_values[i] - y_jupiter_values[i]) / r**3
    ax_sat = -G * M_jupiter * (x_sat_values[i] - x_jupiter_values[i]) / r**3
    ay_sat = -G * M_jupiter * (y_sat_values[i] - y_jupiter_values[i]) / r**3
    
    # Update positions
    x_jupiter_values[i + 1] = x_jupiter_values[i] + vx_jupiter_values[i] * dt + 0.5 * ax_jupiter * dt**2
    y_jupiter_values[i + 1] = y_jupiter_values[i] + vy_jupiter_values[i] * dt + 0.5 * ay_jupiter * dt**2
    x_sat_values[i + 1] = x_sat_values[i] + vx_sat_values[i] * dt + 0.5 * ax_sat * dt**2
    y_sat_values[i + 1] = y_sat_values[i] + vy_sat_values[i] * dt + 0.5 * ay_sat * dt**2
    
    # Calculate new accelerations
    r_new = np.sqrt((x_sat_values[i + 1] - x_jupiter_values[i + 1])**2 + (y_sat_values[i + 1] - y_jupiter_values[i + 1])**2)
    ax_jupiter_new = G * M_satellite * (x_sat_values[i + 1] - x_jupiter_values[i + 1]) / r_new**3
    ay_jupiter_new = G * M_satellite * (y_sat_values[i + 1] - y_jupiter_values[i + 1]) / r_new**3
    ax_sat_new = -G * M_jupiter * (x_sat_values[i + 1] - x_jupiter_values[i + 1]) / r_new**3
    ay_sat_new = -G * M_jupiter * (y_sat_values[i + 1] - y_jupiter_values[i + 1]) / r_new**3
    
    # Update velocities
    vx_jupiter_values[i + 1] = vx_jupiter_values[i] + 0.5 * (ax_jupiter + ax_jupiter_new) * dt
    vy_jupiter_values[i + 1] = vy_jupiter_values[i] + 0.5 * (ay_jupiter + ay_jupiter_new) * dt
    vx_sat_values[i + 1] = vx_sat_values[i] + 0.5 * (ax_sat + ax_sat_new) * dt
    vy_sat_values[i + 1] = vy_sat_values[i] + 0.5 * (ay_sat + ay_sat_new) * dt
    


input_file = open("dat.txt", "w")
for i in range(n_steps):
    input_file.write(f"{t_values[i]} {np.sqrt(vx_sat_values[i]**2+ vy_sat_values[i]**2)} {vx_sat_values[i]} {vy_sat_values[i]} {x_sat_values[i]} {y_sat_values[i]}\n")

input_file.close()

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))

# Set background color to black
fig.patch.set_facecolor('black')
ax1.set_facecolor('black')
ax2.set_facecolor('black')

# Set the aspect ratio and limits
ax1.set_aspect('equal')
scale_lim = 10
ax1.set_xlim(-scale_lim*2.5e9, scale_lim*2.5e9)
ax1.set_ylim(1e9, scale_lim*3e9)

# Plot Jupiter and satellite
jupiter, = ax1.plot([], [], 'o', color='cyan', markersize=10)  # Exaggerated size for visibility
jTrail, = ax1.plot([], [], '-', color='cyan', lw=1, alpha=1)  # Change trail to blue
satellite, = ax1.plot([], [], 'o', color='red')
trail, = ax1.plot([], [], '-', color='red', lw=1, alpha=1)  # Change trail to white

# Customize ax2
ax2.set_xlim(-scale_lim*r0/10, scale_lim*r0/10)
ax2.set_ylim(20000, np.max(np.sqrt(vx_sat_values**2 + vy_sat_values**2)) + 25000)

# Speed plot
speed_line, = ax2.plot([], [], 'r-', label='Speed (m/s)')
ax2.axhline(y=40000, color='dodgerblue', linestyle='--', label='Base Speed') 
ax2.axhline(y=40000 + 13.1e3, color='green', linestyle='--', label='Base Speed + Jupiter\'s Speed')  
ax2.axhline(y=40000 - 13.1e3, color='orange', linestyle='--', label='Base Speed - Jupiter\'s Speed')  
plt.title(label="Orbital Swingby (Collision Case)", fontsize=20, color="white")

# Set white labels and title
ax2.set_xlabel('Position x (m)', color='white')
ax2.set_ylabel('Speed (m/s)', color='white')
ax2.legend(facecolor='black', edgecolor='white', labelcolor='white')

# Change tick colors to white
ax1.tick_params(axis='both', colors='white')
ax2.tick_params(axis='both', colors='white')

# Set the spine colors to white
for spine in ax1.spines.values():
    spine.set_color('white')
for spine in ax2.spines.values():
    spine.set_color('white')

def init():
    jupiter.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    jTrail.set_data([], [])
    speed_line.set_data([], [])
    return jupiter, satellite, trail, speed_line, jTrail



# Load the explosion gif
explosion_gif = Image.open("D:\\Users\\lkhom\\Desktop\\PythonDS\\Orbital-Project\\boom2.gif")


# Explosion display flag
collision_occurred = False
explosion_frame_index = 0


def update(frame):
    global collision_occurred, explosion_frame_index

    frame = (frame * 100) % n_steps
    jupiter.set_data(x_jupiter_values[frame], y_jupiter_values[frame])
    jTrail.set_data(x_jupiter_values[:frame], y_jupiter_values[:frame])
    satellite.set_data(x_sat_values[frame], y_sat_values[frame])
    trail.set_data(x_sat_values[:frame], y_sat_values[:frame])

    # Check for collision
    dist = np.sqrt((x_sat_values[frame] - x_jupiter_values[frame])**2 + (y_sat_values[frame] - y_jupiter_values[frame])**2)
    if dist < R_jupiter and not collision_occurred:
        collision_occurred = True
        explosion_frame_index = 0  # Reset explosion frame index



    # Update speed plot
    speed = np.sqrt(vx_sat_values[frame]**2 + vy_sat_values[frame]**2)
    speed_line.set_data(x_sat_values[:frame], np.sqrt(vx_sat_values[:frame]**2 + vy_sat_values[:frame]**2))

    return jupiter, satellite, trail, speed_line, jTrail

ani = FuncAnimation(fig, update, frames=len(t_values), init_func=init, blit=True, interval=0)

print(f"Starting speed:{np.sqrt(vx_sat_values[0]**2 + vy_sat_values[0]**2)}, Ending speed:{np.sqrt(vx_sat_values[n_steps]**2 + vy_sat_values[n_steps]**2)}")

plt.show()






