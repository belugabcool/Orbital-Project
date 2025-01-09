import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from matplotlib.patches import Circle

# const/cond
G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M_jp = 1898e24  # mass of jp, kg
M_satellite = 1000   # mass of the satellite, kg
R_jp = 71492e3  # radius of jp, m


r0 = R_jp + 1.5e10 
ve0 = 13.1e3  
vs0 = 40e3  
stheta = 32.1  # angle in degrees

init = [
    0, 0,        # x_jp, y_jp
    0, ve0,      # vx_jp, vy_jp
    -r0, -r0*0.25,  # x_sat, y_sat
    vs0*np.cos(np.radians(stheta)), vs0*np.sin(np.radians(stheta))
]

# time para
t0 = 0
tf = (5*r0)/vs0
n_steps = 50000
dt = (tf - t0) / n_steps

t = np.linspace(t0, tf, n_steps + 1)
x_jp = np.zeros(n_steps + 1)
y_jp = np.zeros(n_steps + 1)
vx_jp = np.zeros(n_steps + 1)
vy_jp = np.zeros(n_steps + 1)
x_sat = np.zeros(n_steps + 1)
y_sat = np.zeros(n_steps + 1)
vx_sat = np.zeros(n_steps + 1)
vy_sat = np.zeros(n_steps + 1)
fig.patch.set_facecolor('black')
ax1.set_facecolor('black')
ax2.set_facecolor('black')

x_jp[0], y_jp[0] = init[0], init[1]
vx_jp[0], vy_jp[0] = init[2], init[3]
x_sat[0], y_sat[0] = init[4], init[5]
vx_sat[0], vy_sat[0] = init[6], init[7]


collision_check = False
collision_frame = None  

# velocity verlet (newton grav)
# ax = g * m * (delta x) / r^3 
# ay = g * m * (delta y) / r^3 
# x = x + vx * timestep + 0.5 * ax * time step^2 
# y = y + yx * timestep + 0.5 * ay * time step^2 
# delta v = 0.5 * (a[t] + a[t+1]) * time step 
for i in range(n_steps):
    r = np.sqrt(
        (x_sat[i] - x_jp[i])**2 
        + (y_sat[i] - y_jp[i])**2
    )
    
    # collision check
    if r < R_jp and not collision_check:
        collision_check = True
        collision_frame = i

    if collision_check:

        x_jp[i + 1] = x_jp[i]
        y_jp[i + 1] = y_jp[i] 
        vx_jp[i + 1] = 0
        vy_jp[i + 1] = 0

        x_sat[i + 1] = x_sat[i]
        y_sat[i + 1] = y_sat[i]
        vx_sat[i + 1] = 0
        vy_sat[i + 1] = 0
        continue 
    

    ax_jp = (G * M_satellite * (x_sat[i] - x_jp[i]) / r**3)
    ay_jp = (G * M_satellite * (y_sat[i] - y_jp[i]) / r**3)
    ax_sat = (-G * M_jp * (x_sat[i] - x_jp[i]) / r**3)
    ay_sat = (-G * M_jp * (y_sat[i] - y_jp[i]) / r**3)
    

    x_jp[i + 1] = (x_jp[i] + vx_jp[i] * dt + 0.5 * ax_jp * dt**2)
    y_jp[i + 1] = (y_jp[i] + vy_jp[i] * dt + 0.5 * ay_jp * dt**2)
    x_sat[i + 1] = (x_sat[i] + vx_sat[i] * dt + 0.5 * ax_sat * dt**2)

    y_sat[i + 1] = (y_sat[i] + vy_sat[i] * dt + 0.5 * ay_sat * dt**2)
    

    r_new = np.sqrt((x_sat[i + 1] - x_jp[i + 1])**2+ (y_sat[i + 1] - y_jp[i + 1])**2)
    ax_jp_new = (G * M_satellite * (x_sat[i + 1] - x_jp[i + 1]) / r_new**3)
    ay_jp_new = (G * M_satellite * (y_sat[i + 1] - y_jp[i + 1]) / r_new**3)
    ax_sat_new = (-G * M_jp * (x_sat[i + 1] - x_jp[i + 1]) / r_new**3)
    ay_sat_new = (-G * M_jp * (y_sat[i + 1] - y_jp[i + 1]) / r_new**3)
    

    vx_jp[i + 1] = (vx_jp[i] + 0.5 * (ax_jp + ax_jp_new) * dt)
    vy_jp[i + 1] = (vy_jp[i] + 0.5 * (ay_jp + ay_jp_new) * dt)
    vx_sat[i + 1] = (vx_sat[i] + 0.5 * (ax_sat + ax_sat_new) * dt)
    vy_sat[i + 1] = (vy_sat[i] + 0.5 * (ay_sat + ay_sat_new) * dt)

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 12))




ax1.set_aspect('equal')
scale_lim = 10
ax1.set_xlim(-scale_lim*2.5e9, scale_lim*2.5e9)
ax1.set_ylim(1e9, scale_lim*3e9)

jp, = ax1.plot([], [], 'o', color='cyan', markersize=10, label='jp')  
jTrail, = ax1.plot([], [], '-', color='cyan', lw=1, alpha=1)
satellite, = ax1.plot([], [], 'o', color='red', label='Satellite')
trail, = ax1.plot([], [], '-', color='red', lw=1, alpha=1)


explosion_red = Circle((0, 0), 0, color='red', alpha=1, zorder=10)
explosion_orange = Circle((0, 0), 0, color='orange', alpha=0.6, zorder=10)
explosion_yellow = Circle((0, 0), 0, color='yellow', alpha=0.4, zorder=10)


ax1.add_patch(explosion_red)
ax1.add_patch(explosion_orange)
ax1.add_patch(explosion_yellow)
ax1.add_patch(explosion_salmon)
ax1.add_patch(explosion_firebrick)
ax1.add_patch(explosion_ired)

ax2.set_xlim(-scale_lim*r0/10, scale_lim*r0/10)
ax2.set_ylim(
    20000, 
    np.max(np.sqrt(vx_sat**2 + vy_sat**2)) + 25000
)
speed_line, = ax2.plot([], [], 'r-', label='Speed (m/s)')
ax2.axhline(y=40000, color='dodgerblue', linestyle='--', label='Base Speed')
ax2.axhline(y=40000 + 13.1e3, color='green', linestyle='--', 
            label='Base Speed + jp\'s Speed')
ax2.axhline(y=40000 - 13.1e3, color='orange', linestyle='--', 
            label='Base Speed - jp\'s Speed')

ax2.set_xlabel('Position x (m)', color='white')
ax2.set_ylabel('speed (m/s)', color='white')
ax2.legend(facecolor='black', edgecolor='white', labelcolor='white')
plt.title(label="Orbital Sswingby (Collision Case)", fontsize=20, color="white")

ax1.tick_params(axis='both', colors='white')
ax2.tick_params(axis='both', colors='white')

for spine in ax1.spines.values():
    spine.set_color('white')
for spine in ax2.spines.values():
    spine.set_color('white')

collision_occurred = False  
explosion_start_frame = None


collision_occurred = False
explosion_start_frame = None
prev_frame = -1  
explosion_salmon = Circle((0, 0), 0, color='ghostwhite', alpha=0.3, zorder=10)
explosion_firebrick = Circle((0, 0), 0, color='darkorange', alpha=0.5, zorder=10)
explosion_ired = Circle((0, 0), 0, color='indianred', alpha=0.6, zorder=10)
def init():

    jp.set_data([], [])
    jTrail.set_data([], [])
    satellite.set_data([], [])
    trail.set_data([], [])
    speed_line.set_data([], [])
    satellite.set_visible(True)

    explosion_red.set_radius(0)
    explosion_orange.set_radius(0)
    explosion_yellow.set_radius(0)
    explosion_salmon.set_radius(0)
    explosion_firebrick.set_radius(0)
    explosion_ired.set_radius(0)
    return (
        jp, jTrail, satellite, trail,speed_line, explosion_red, explosion_orange, explosion_yellow, explosion_salmon,explosion_firebrick, explosion_ired)

explosioncount = 0 

def update(frame):

    global collision_occurred, explosion_start_frame, prev_frame, explosioncount

    

    frame = (frame * 100) % n_steps


    jp.set_data(x_jp[frame], y_jp[frame])
    jTrail.set_data(x_jp[:frame], y_jp[:frame])
    satellite.set_data(x_sat[frame], y_sat[frame])
    trail.set_data(x_sat[:frame], y_sat[:frame])

    #  collision
    dist = np.sqrt(
        (x_sat[frame] - x_jp[frame])**2 
        + (y_sat[frame] - y_jp[frame])**2
    )
    if dist < R_jp and not collision_occurred:
        collision_occurred = True
        explosion_start_frame = frame
        print("hi")

        explosioncount += 1
        if explosioncount != 1:
            print(explosioncount)
            
        

    if collision_occurred:
        

        satellite.set_visible(True)

        cx = x_jp[frame]
        cy = y_jp[frame]

        frames_since_collision = frame - explosion_start_frame
        if frames_since_collision < 0:
            frames_since_collision = 0

        if frames_since_collision >= 100:
            satellite.set_visible(False)
            collision_occured = False


        explosion_red.set_center((cx, cy))
        explosion_orange.set_center((cx, cy))
        explosion_yellow.set_center((cx, cy))
        explosion_salmon.set_center((cx, cy))
        explosion_firebrick.set_center((cx, cy))
        explosion_ired.set_center((cx, cy))


        explosion_red.set_radius(100*0.0004/2 * frames_since_collision * R_jp)
        explosion_orange.set_radius(100*0.0007/2 * frames_since_collision  * R_jp)
        explosion_yellow.set_radius(100*0.0010/2 *  abs(16*frames_since_collision - (frames_since_collision/100)**1.74 * 50)/3 * R_jp)
        explosion_salmon.set_radius(100*0.0007/2 * R_jp * abs(16*frames_since_collision - (frames_since_collision/100)**1.73 * 50)/3)
        explosion_firebrick.set_radius(100*0.0010/2   * R_jp * frames_since_collision * abs(16*frames_since_collision - (frames_since_collision/100)**1.7 * 50)/3)
        explosion_ired.set_radius(100*0.0010/2   * R_jp * frames_since_collision * abs(16*frames_since_collision - (frames_since_collision/100)**1.7 * 20)/3)
        print("a + " + str(frames_since_collision))
        print(16*frames_since_collision) 
        print( (frames_since_collision/100)**1.8 * 100)
        print(16*frames_since_collision - (frames_since_collision/100)**1.8 * 100)



        fade_factor = 1 - frames_since_collision / 30000
        fade_factor = max(fade_factor, 0)
        explosion_red.set_alpha(0.8 * fade_factor)
        explosion_orange.set_alpha(0.6 * fade_factor)
        explosion_yellow.set_alpha(0.4 * fade_factor)
        explosion_salmon.set_alpha(0.3 * fade_factor)
        explosion_firebrick.set_alpha(0.5 * fade_factor)




    speed_line.set_data(
        x_sat[:frame], 
        np.sqrt(vx_sat[:frame]**2 + vy_sat[:frame]**2)
    )

    return (jp, jTrail, satellite, trail,
        speed_line, explosion_red, explosion_orange, explosion_yellow, explosion_salmon, explosion_firebrick, explosion_ired)

ani = FuncAnimation(
    fig, 
    update, 
    frames=len(t), 
    init_func=init, 
    blit=True, 
    interval=0
)

print(
    f"Starting speed: {np.sqrt(vx_sat[0]**2 + vy_sat[0]**2):.2f} m/s, "
    f"Ending speed:   {np.sqrt(vx_sat[n_steps]**2 + vy_sat[n_steps]**2):.2f} m/s"
)


plt.show()
