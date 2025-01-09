import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

# Constants
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_jupiter = 1898e24  # mass of Jupiter, kg
M_satellite = 1000   # mass of the satellite, kg
R_jupiter = 71492e3  # radius of Jupiter, m

# Initial conditions
r0 = R_jupiter + 1.5e10 
ve0 = 13.1e3  
vs0 = 40e3  
stheta = 32.1  # angle in degrees

initial_state = [
    0, 0,        # x_jupiter, y_jupiter
    0, ve0,      # vx_jupiter, vy_jupiter
    -r0, -r0*0.25,  # x_sat, y_sat
    vs0*np.cos(np.radians(stheta)), vs0*np.sin(np.radians(stheta))
]

# Time parameters
t0 = 0
tf = (5*r0)/vs0
n_steps = 50000
dt = (tf - t0) / n_steps

# Arrays to store the solution
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

# Velocity Verlet to solve equations of motion
collision_check = False
collision_frame = None  # Will store the frame index at collision

for i in range(n_steps):
    r = np.sqrt(
        (x_sat_values[i] - x_jupiter_values[i])**2 
        + (y_sat_values[i] - y_jupiter_values[i])**2
    )
    
    # Check collision
    if r < R_jupiter and not collision_check:
        collision_check = True
        collision_frame = i

    if collision_check:
        # Once collision has happened, freeze everything
        x_jupiter_values[i + 1] = x_jupiter_values[i]
        y_jupiter_values[i + 1] = y_jupiter_values[i] 
        vx_jupiter_values[i + 1] = 0
        vy_jupiter_values[i + 1] = 0

        x_sat_values[i + 1] = x_sat_values[i]
        y_sat_values[i + 1] = y_sat_values[i]
        vx_sat_values[i + 1] = 0
        vy_sat_values[i + 1] = 0
        continue 
    
    # Accelerations
    ax_jupiter = (
        G * M_satellite 
        * (x_sat_values[i] - x_jupiter_values[i]) / r**3
    )
    ay_jupiter = (
        G * M_satellite 
        * (y_sat_values[i] - y_jupiter_values[i]) / r**3
    )
    ax_sat = (
        -G * M_jupiter 
        * (x_sat_values[i] - x_jupiter_values[i]) / r**3
    )
    ay_sat = (
        -G * M_jupiter 
        * (y_sat_values[i] - y_jupiter_values[i]) / r**3
    )
    
    # Update positions
    x_jupiter_values[i + 1] = (
        x_jupiter_values[i] 
        + vx_jupiter_values[i] * dt 
        + 0.5 * ax_jupiter * dt**2
    )
    y_jupiter_values[i + 1] = (
        y_jupiter_values[i] 
        + vy_jupiter_values[i] * dt 
        + 0.5 * ay_jupiter * dt**2
    )
    x_sat_values[i + 1] = (
        x_sat_values[i] 
        + vx_sat_values[i] * dt 
        + 0.5 * ax_sat * dt**2
    )
    y_sat_values[i + 1] = (
        y_sat_values[i] 
        + vy_sat_values[i] * dt 
        + 0.5 * ay_sat * dt**2
    )
    
    # Calculate new accelerations
    r_new = np.sqrt(
        (x_sat_values[i + 1] - x_jupiter_values[i + 1])**2
        + (y_sat_values[i + 1] - y_jupiter_values[i + 1])**2
    )
    ax_jupiter_new = (
        G * M_satellite 
        * (x_sat_values[i + 1] - x_jupiter_values[i + 1]) / r_new**3
    )
    ay_jupiter_new = (
        G * M_satellite 
        * (y_sat_values[i + 1] - y_jupiter_values[i + 1]) / r_new**3
    )
    ax_sat_new = (
        -G * M_jupiter 
        * (x_sat_values[i + 1] - x_jupiter_values[i + 1]) / r_new**3
    )
    ay_sat_new = (
        -G * M_jupiter 
        * (y_sat_values[i + 1] - y_jupiter_values[i + 1]) / r_new**3
    )
    
    # Update velocities
    vx_jupiter_values[i + 1] = (
        vx_jupiter_values[i] 
        + 0.5 * (ax_jupiter + ax_jupiter_new) * dt
    )
    vy_jupiter_values[i + 1] = (
        vy_jupiter_values[i] 
        + 0.5 * (ay_jupiter + ay_jupiter_new) * dt
    )
    vx_sat_values[i + 1] = (
        vx_sat_values[i] 
        + 0.5 * (ax_sat + ax_sat_new) * dt
    )
    vy_sat_values[i + 1] = (
        vy_sat_values[i] 
        + 0.5 * (ay_sat + ay_sat_new) * dt
    )

# --- PLOTTING ---
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))

# Set background color to black
fig.patch.set_facecolor('black')
ax1.set_facecolor('black')
ax2.set_facecolor('black')

# Set aspect ratio and limits
ax1.set_aspect('equal')
scale_lim = 10
ax1.set_xlim(-scale_lim*2.5e9, scale_lim*2.5e9)
ax1.set_ylim(1e9, scale_lim*3e9)

# Plot Jupiter and satellite
jupiter, = ax1.plot([], [], 'o', color='cyan', markersize=10, label='Jupiter')  
jTrail, = ax1.plot([], [], '-', color='cyan', lw=1, alpha=1)
satellite, = ax1.plot([], [], 'o', color='red', label='Satellite')
trail, = ax1.plot([], [], '-', color='red', lw=1, alpha=1)

# Explosion circles (initially radius=0)
explosion_red = Circle((0, 0), 0, color='red', alpha=1, zorder=10)
explosion_orange = Circle((0, 0), 0, color='orange', alpha=0.6, zorder=10)
explosion_yellow = Circle((0, 0), 0, color='yellow', alpha=0.4, zorder=10)
explosion_salmon = Circle((0, 0), 0, color='ghostwhite', alpha=0.3, zorder=10)
explosion_firebrick = Circle((0, 0), 0, color='darkorange', alpha=0.5, zorder=10)
explosion_ired = Circle((0, 0), 0, color='indianred', alpha=0.6, zorder=10)

ax1.add_patch(explosion_red)
ax1.add_patch(explosion_orange)
ax1.add_patch(explosion_yellow)
ax1.add_patch(explosion_salmon)
ax1.add_patch(explosion_firebrick)
ax1.add_patch(explosion_ired)
# Speed plot
ax2.set_xlim(-scale_lim*r0/10, scale_lim*r0/10)
ax2.set_ylim(
    20000, 
    np.max(np.sqrt(vx_sat_values**2 + vy_sat_values**2)) + 25000
)
speed_line, = ax2.plot([], [], 'r-', label='Speed (m/s)')
ax2.axhline(y=40000, color='dodgerblue', linestyle='--', label='Base Speed')
ax2.axhline(y=40000 + 13.1e3, color='green', linestyle='--', 
            label='Base Speed + Jupiter\'s Speed')
ax2.axhline(y=40000 - 13.1e3, color='orange', linestyle='--', 
            label='Base Speed - Jupiter\'s Speed')

ax2.set_xlabel('Position x (m)', color='white')
ax2.set_ylabel('Speed (m/s)', color='white')
ax2.legend(facecolor='black', edgecolor='white', labelcolor='white')
plt.title(label="Orbital Swingby (Collision Case)", fontsize=20, color="white")
# Change tick colors to white
ax1.tick_params(axis='both', colors='white')
ax2.tick_params(axis='both', colors='white')
# Spine colors to white
for spine in ax1.spines.values():
    spine.set_color('white')
for spine in ax2.spines.values():
    spine.set_color('white')

collision_occurred = False  # We'll also track it here
explosion_start_frame = None

# Globals to track collision and looping
collision_occurred = False
explosion_start_frame = None
prev_frame = -1  # store previous frame so we know if we've looped

def init():
    """Initialize animation frames."""
    jupiter.set_data([], [])
    jTrail.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    speed_line.set_data([], [])
    satellite.set_visible(True)
    # Hide/zero the explosion circles
    explosion_red.set_radius(0)
    explosion_orange.set_radius(0)
    explosion_yellow.set_radius(0)
    explosion_salmon.set_radius(0)
    explosion_firebrick.set_radius(0)
    explosion_ired.set_radius(0)
    return (
        jupiter, jTrail, satellite, trail,
        speed_line, explosion_red, explosion_orange, explosion_yellow, explosion_salmon,explosion_firebrick, explosion_ired
    )

explosioncount = 0 

def update(frame):
    """Update function for FuncAnimation."""
    global collision_occurred, explosion_start_frame, prev_frame, explosioncount

    
    # Step in bigger increments to shorten animation time
    frame = (frame * 100) % n_steps

    # Update Jupiter & Satellite paths
    jupiter.set_data(x_jupiter_values[frame], y_jupiter_values[frame])
    jTrail.set_data(x_jupiter_values[:frame], y_jupiter_values[:frame])
    satellite.set_data(x_sat_values[frame], y_sat_values[frame])
    trail.set_data(x_sat_values[:frame], y_sat_values[:frame])

    # Check collision
    dist = np.sqrt(
        (x_sat_values[frame] - x_jupiter_values[frame])**2 
        + (y_sat_values[frame] - y_jupiter_values[frame])**2
    )
    if dist < R_jupiter and not collision_occurred:
        collision_occurred = True
        explosion_start_frame = frame
        print("hi")

        explosioncount += 1
        if explosioncount != 1:
            print(explosioncount)
            
        

    if collision_occurred:
        

        # Hide the satellite once collision occurs
        satellite.set_visible(True)
        # Position of the explosion = Jupiter's center
        cx = x_jupiter_values[frame]
        cy = y_jupiter_values[frame]

        # How many frames since collision started
        frames_since_collision = frame - explosion_start_frame
        if frames_since_collision < 0:
            frames_since_collision = 0

        if frames_since_collision >= 100:
            satellite.set_visible(False)
            collision_occured = False

        # Expand circles with time (tweak scale factors as desired)
        explosion_red.set_center((cx, cy))
        explosion_orange.set_center((cx, cy))
        explosion_yellow.set_center((cx, cy))
        explosion_salmon.set_center((cx, cy))
        explosion_firebrick.set_center((cx, cy))
        explosion_ired.set_center((cx, cy))


        explosion_red.set_radius(100*0.0004/2 * frames_since_collision * R_jupiter)
        explosion_orange.set_radius(100*0.0007/2 * frames_since_collision  * R_jupiter)
        explosion_yellow.set_radius(100*0.0010/2 *  abs(16*frames_since_collision - (frames_since_collision/100)**1.74 * 50)/3 * R_jupiter)
        explosion_salmon.set_radius(100*0.0007/2 * R_jupiter * abs(16*frames_since_collision - (frames_since_collision/100)**1.73 * 50)/3)
        explosion_firebrick.set_radius(100*0.0010/2   * R_jupiter * frames_since_collision * abs(16*frames_since_collision - (frames_since_collision/100)**1.7 * 50)/3)
        explosion_ired.set_radius(100*0.0010/2   * R_jupiter * frames_since_collision * abs(16*frames_since_collision - (frames_since_collision/100)**1.7 * 20)/3)
        print("a + " + str(frames_since_collision))
        print(16*frames_since_collision) 
        print( (frames_since_collision/100)**1.8 * 100)
        print(16*frames_since_collision - (frames_since_collision/100)**1.8 * 100)


        # Optionally fade them out over time
        # e.g., alpha decreases from 1.0 to 0.0 in 300 frames
        fade_factor = 1 - frames_since_collision / 30000
        fade_factor = max(fade_factor, 0)
        explosion_red.set_alpha(0.8 * fade_factor)
        explosion_orange.set_alpha(0.6 * fade_factor)
        explosion_yellow.set_alpha(0.4 * fade_factor)
        explosion_salmon.set_alpha(0.3 * fade_factor)
        explosion_firebrick.set_alpha(0.5 * fade_factor)




    # Update speed plot
    speed_line.set_data(
        x_sat_values[:frame], 
        np.sqrt(vx_sat_values[:frame]**2 + vy_sat_values[:frame]**2)
    )

    return (
        jupiter, jTrail, satellite, trail,
        speed_line, explosion_red, explosion_orange, explosion_yellow, explosion_salmon, explosion_firebrick, explosion_ired
    )

ani = FuncAnimation(
    fig, 
    update, 
    frames=len(t_values), 
    init_func=init, 
    blit=True, 
    interval=0
)

print(
    f"Starting speed: {np.sqrt(vx_sat_values[0]**2 + vy_sat_values[0]**2):.2f} m/s, "
    f"Ending speed:   {np.sqrt(vx_sat_values[n_steps]**2 + vy_sat_values[n_steps]**2):.2f} m/s"
)








plt.show()
