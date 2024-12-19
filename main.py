import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.special as sps
from matplotlib import animation
import ipywidgets as widg
from IPython.display import HTML
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk
from functions_and_classes import *

# Create the grids for every option

# Rectangle membrane dimensions
Lx = 1
Ly = 1

# Circular membrane dimensions
radius = 1

# Grid for rectangle membrane with boundary condition 0_onS
x0 = np.linspace(0, Lx, 250)
y0 = np.linspace(0, Ly, 250)
xv0, yv0 = np.meshgrid(x0, y0)

# Grid for rectangle membrane with boundary condition Max_onS
x = np.linspace(-Lx/2, Ly/2, 250)
y = np.linspace(-Lx/2, Ly/2, 250)
xv, yv = np.meshgrid(x, y) 


# Grid for circular membrane (symmetric and general) with boundary condition 0_onS
rj = np.linspace(0, radius, 250) # radius=1
thetaj= np.linspace(0, 2*np.pi, 250)
rjv, thetajv = np.meshgrid(rj, thetaj)



# Initialize figure and axes based on membrane shape and animation type (waves or particles, but for now only particles)
class Membrane():
    """_summary_
    """
    def __init__(self, mtype='rect', boundary='0_onS', atype='particles') -> None:
        self.type = mtype
        self.boundary = boundary
        self.animation_type = atype
        self.resolution = 1 # for slider intervals
        self.a = 0 # for uniform distribution
        self.b = 1
        

        # Setting initial amplitude in the membrane
        if self.type == 'rect':
            self.fig, self.ax = plt.subplots(figsize=(8,6))
            if self.boundary == '0_onS':
                self.ax.set_xlim(0, Lx)
                self.ax.set_ylim(0, Ly)
                self.x = xv0
                self.y = yv0
                self.amplitude = wave_amp(x_vals=self.x, y_vals=self.y, n=1, m=3, L_x=Lx, L_y=Ly, boundary=self.boundary, mtype=self.type, abso='no') # for waves animation
                self.a = 0 # for uniform distribution
                self.b = Lx
                self.omega = rect_omega(n=1, m=3, L_x=Lx, L_y=Ly, wave_speed=0.75)
            else:
                self.ax.set_xlim(-Lx/2, Lx/2)
                self.ax.set_ylim(-Ly/2, Ly/2)
                self.x = xv
                self.y = yv
                self.amplitude = wave_amp(x_vals=self.x, y_vals=self.y, n=1, m=3, L_x=Lx, L_y=Ly, boundary=self.boundary, mtype=self.type, abso='no') # for waves animation
                self.resolution = 2 # for slider intervals to get only odd numbers
                self.a = -Lx/2 # for uniform distribution
                self.b = Lx/2
                self.omega = rect_omega(n=1, m=3, L_x=Lx, L_y=Ly, wave_speed=0.75)


        elif 'circ' in self.type:
            self.fig, self.ax = plt.subplots(subplot_kw={'projection':'polar'}, figsize=(8,6))
            self.ax.set_ylim(0, radius + 0.05) # radius value + 0.05 for better visualization
            self.x = rjv
            self.y = thetajv
            self.boundary = '0_onS'
            self.amplitude = wave_amp(x_vals=self.x, y_vals=self.y, n=1, m=3, boundary='0_onS', mtype=self.type, abso='no') # for waves animation
            self.a = 0 # for uniform distribution
            self.b = 2*np.pi
            self.omega = circ_omega(n=1, m=3, radius=radius, wave_speed=0.75)
        else:
            print("Please select one of the following options: 'rect', 'circ_sym' or 'circ_gen'.")
            self.amplitude = 0

        # Set ticks to be invisible
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)

        # Setting plot based on animation type
        if self.animation_type == 'particles':
            self.plot_for_particles()
        elif self.animation_type == 'waves':
            self.plot_for_waves()


    def plot_for_particles(self) -> None:
        self.dot, = plt.plot([],[], 'o', ms=2, color='#cbbd93') # sand color #cbbd93
        self.ax.set_facecolor('dimgrey')    
    

    def plot_for_waves(self) -> None:
        # self.fig.set_size_inches(10,8)
        self.ax.axis('off')
        self.update_wave_ax(self.amplitude)
        im = self.ax.imshow(self.amplitude)
        self.fig.colorbar(mappable=im, ax=self.ax)
    

    def update_wave_ax(self, amplitude) -> None:
        self.ax.clear()
        if 'circ' in self.type:
            self.ax.pcolormesh(self.y, self.x, amplitude)
        else:
            self.ax.pcolormesh(self.x, self.y, amplitude)

        





def main():
    # Tkinter GUI setup
    root = tk.Tk()
    root.title("Interactive Chladni Figures")
    root.minsize(1000, 800)

    def set_global_variables():
        global waves_or_particles, rect_or_circ, zero_or_max_onS, wave_speed_in_membrane, absolute
        waves_or_particles = 'particles'
        rect_or_circ = 'rect'
        zero_or_max_onS = '0_onS'

        wave_speed_in_membrane = 0.75
        absolute = 'no'

    set_global_variables()

    # Update function for the animation,
    def animate(i):
        global ensemble, waves_or_particles, wave_speed_in_membrane, absolute
        n_mode = slider_n.get()
        m_mode = slider_m.get()
        if 'membrane' in globals():
            if waves_or_particles == 'particles':
                if i%5 == 0:
                    ensemble.step(L_x=Lx, L_y=Ly, n=n_mode, m=m_mode, boundary=membrane.boundary, mtype=membrane.type, abso=absolute)

                points = ensemble.prev_points + (i%5)/5 * (ensemble.points - ensemble.prev_points)
                if 'circ' in membrane.type:
                    membrane.dot.set_data(points[1], points[0])
                else:
                    membrane.dot.set_data(*points)
            elif waves_or_particles == 'waves':
                if membrane.type == 'rect':
                    omega = rect_omega(n=n_mode, m=m_mode, L_x=Lx, L_y=Ly, wave_speed=wave_speed_in_membrane)
                elif 'circ' in membrane.type:
                    omega = circ_omega(n=n_mode, m=m_mode, radius=radius, wave_speed=wave_speed_in_membrane)
                amplitude = wave_amp(x_vals=membrane.x, y_vals=membrane.y, n=n_mode, m=m_mode, L_x=Lx, L_y=Ly, 
                                     boundary=membrane.boundary, mtype=membrane.type, abso=absolute) * time_evolution(omega=omega, t=i / 10)
                membrane.update_wave_ax(amplitude=amplitude)
                

    # Create function for reset particle animation
    def create_particles():
        global ensemble
        if 'membrane' in globals():
            ensemble = Particles(amplitude=wave_amp, a=membrane.a, b=membrane.b, ttype=membrane.type, num_points=10000, delta=0.1)
    

    # Initialization function for start animation
    def start_animation():
        global ani
        ani = animation.FuncAnimation(membrane.fig, animate, frames=500, interval=30, repeat=True)
        canvas.draw_idle()  # ensure the canvas updates


    # Stop function for stop animation
    def stop_animation():
        global ani
        if 'ani' in globals():
            ani.event_source.stop()  # stop the animation
    

    # Radiobuttons functions
    # Check the values for animation type, shape and boundary. It should be the same for the membrane attributes and global variables.
    def test_changes(*args):
        print(f'\nanimation radiobutton: {waves_or_particles}, membrane animation: {membrane.animation_type}')
        print(f'shape radiobutton: {rect_or_circ}, membrane shape: {membrane.type}')
        print(f'boundary radiobutton: {zero_or_max_onS}, membrane boundary: {membrane.boundary}')
        print(f'Is the wave amplitude on absolute value? {absolute}')
        print(f'The wave speed for wave animation is: {wave_speed_in_membrane}')


    # Select animation function
    def select_animation():
        # global waves_or_particles
        global waves_or_particles, rect_or_circ, zero_or_max_onS, absolute
        if animation_var.get() == 1:
            waves_or_particles = 'particles'
            if reset_button.cget("state") == "disabled": # activate the reset particles button if disabled
                reset_button.config(state="normal")
            if speed_button.cget("state") == "normal": # disable the wave speed apply button if activated
                speed_button.config(state="disabled")
        elif animation_var.get() == 2:
            waves_or_particles = 'waves'
            if reset_button.cget("state") == "normal": # disable the reset particles button for waves animation
                reset_button.config(state="disabled")
            if speed_button.cget("state") == "disabled": # activate the wave speed apply button if disabled
                speed_button.config(state="normal")

        initiate_membrane(rect_or_circ, zero_or_max_onS, waves_or_particles)
        test_changes(rect_or_circ, zero_or_max_onS, waves_or_particles)
        # return waves_or_particles
    

    # Select membrane shape
    def select_shape():
        # global rect_or_circ
        global waves_or_particles, rect_or_circ, zero_or_max_onS
        if shape_var.get() == 1:
            rect_or_circ = 'rect'
        elif shape_var.get() == 2:
            rect_or_circ = 'circ_sym'
            boundary_var.set(1)
            select_boundary()
        elif shape_var.get() == 3:
            rect_or_circ = 'circ_gen'
            boundary_var.set(1)
            select_boundary()
        slider_n.set(1)
        slider_m.set(3)
        initiate_membrane(rect_or_circ, zero_or_max_onS, waves_or_particles)
        test_changes(rect_or_circ, zero_or_max_onS, waves_or_particles)
        # return rect_or_circ


    # Select boundary condition
    def select_boundary():
        # global zero_or_max_onS
        global waves_or_particles, rect_or_circ, zero_or_max_onS
        if boundary_var.get() == 1:
            zero_or_max_onS = '0_onS'
            slider_n.set(1)
            slider_m.set(3)
            slider_n.config(resolution=1)
            slider_m.config(resolution=1)
        elif boundary_var.get() == 2:
            zero_or_max_onS = 'Max_onS'
            shape_var.set(1)
            select_shape()
            slider_n.set(1)
            slider_m.set(3)
            slider_n.config(resolution=2)
            slider_m.config(resolution=2)
        initiate_membrane(rect_or_circ, zero_or_max_onS, waves_or_particles)
        test_changes(rect_or_circ, zero_or_max_onS, waves_or_particles)
        # return zero_or_max_onS


    # Select absolute value
    def select_absolute():
        global absolute
        if absolute_var.get() == 1:
            absolute = 'no'
        elif absolute_var.get() == 2:
            absolute = 'yes'
        test_changes()


    def change_speed():
        global wave_speed_in_membrane
        wave_speed_in_membrane = float(speed_entry.get())

    # Initialize the membrane
    def initiate_membrane(rect_or_circ, zero_or_max_onS, waves_or_particles):
        global membrane
        stop_animation() 
        membrane = Membrane(mtype=rect_or_circ, boundary=zero_or_max_onS, atype=waves_or_particles) # for mtype: 'rect', 'circ_sym' or 'circ_gen'. for atype: 'particles' or 'waves'.
        create_particles() # initialize the particles
        if 'canvas' in globals(): 
            canvas.get_tk_widget().pack_forget() 
            initiate_canves()
        

    # Initialize canvas
    def initiate_canves():
        # Main window
        global canvas
        if 'membrane' in globals():
            canvas = FigureCanvasTkAgg(membrane.fig, master=frame)
            canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)
        

    # Control panels
    frame = tk.Frame(root)
    frame.pack(fill=tk.BOTH, expand=1)

    right_frame = tk.Frame(frame)
    right_frame.pack(side=tk.RIGHT, fill=tk.Y)
    
    bottom_frame = tk.Frame(root)
    bottom_frame.pack(side=tk.BOTTOM, fill=tk.X)

    # Initialize objects
    initiate_membrane(rect_or_circ, zero_or_max_onS, waves_or_particles) # initialize the membrane
    initiate_canves() # initialize the canvas


    # Sliders for n and m
    label_n = tk.Label(bottom_frame, text='n : ')
    label_n.pack(side=tk.LEFT, padx=(100, 0))
    slider_n = tk.Scale(bottom_frame, from_=1, to=15, resolution=membrane.resolution, orient="horizontal")
    slider_n.set(1)
    slider_n.pack(side=tk.LEFT, padx=(5, 10))

    label_m = tk.Label(bottom_frame, text='m : ')
    label_m.pack(side=tk.LEFT, padx=(10, 0))
    slider_m = tk.Scale(bottom_frame, from_=1, to=15, resolution=membrane.resolution, orient="horizontal")
    slider_m.set(3)
    slider_m.pack(side=tk.LEFT, padx=(5, 10))

    # Buttons for start, stop and reset animation
    start_button = tk.Button(bottom_frame, text="Start Animation", command=start_animation)
    start_button.pack(side=tk.LEFT, padx=10)

    stop_button = tk.Button(bottom_frame, text="Stop Animation", command=stop_animation)
    stop_button.pack(side=tk.LEFT, padx=10)

    reset_button = tk.Button(bottom_frame, text="Reset Particles", command=create_particles) # , command=create_particles
    reset_button.pack(side=tk.LEFT, padx=10)


    # Create radiobuttons for mulitple options for particle/plane, rectangular/circular, 0_onS/Max_onS

    # Radiobuttons for particle/plane
    animation_var = tk.IntVar(master=right_frame, value=1)

    animation_label = tk.Label(master=right_frame, text='Animation:')
    animation_label.pack(fill=tk.X, padx=10, pady=(30, 0))

    particle_option = tk.Radiobutton(master=right_frame, text='Particles', variable=animation_var, value=1, indicator=0, background="light blue", command=select_animation)
    particle_option.pack(fill=tk.X, padx=10, pady=(10, 5))

    waves_option = tk.Radiobutton(master=right_frame, text='Waves', variable=animation_var, value=2, indicator=0, background="light blue", command=select_animation)
    waves_option.pack(fill=tk.X, padx=10, pady=(10))


    # Radiobuttons for rectangular/circular(sym or gen) membrane shape
    shape_var = tk.IntVar(master=right_frame, value=1)

    shape_label = tk.Label(master=right_frame, text='Shape:')
    shape_label.pack(fill=tk.X, padx=10, pady=(30, 0))

    rectangular_option = tk.Radiobutton(master=right_frame, text='rect', variable=shape_var, value=1, indicator=0, background="light blue", command=select_shape)
    rectangular_option.pack(fill=tk.X, padx=10, pady=(10, 5))

    circular_sym_option = tk.Radiobutton(master=right_frame, text='circ_sym', variable=shape_var, value=2, indicator=0, background="light blue", command=select_shape)
    circular_sym_option.pack(fill=tk.X, padx=10, pady=(10, 5))

    circular_gen_option = tk.Radiobutton(master=right_frame, text='circ_gen', variable=shape_var, value=3, indicator=0, background="light blue", command=select_shape)
    circular_gen_option.pack(fill=tk.X, padx=10, pady=(10))


    # Radiobuttons for 0_onS/Max_onS boundary condition
    boundary_var = tk.IntVar(master=right_frame, value=1)

    boundary_label = tk.Label(master=right_frame, text='Boundary:')
    boundary_label.pack(fill=tk.X, padx=10, pady=(30, 0))

    zero_onS_option = tk.Radiobutton(master=right_frame, text='0_onS', variable=boundary_var, value=1, indicator=0, background="light blue", command=select_boundary)
    zero_onS_option.pack(fill=tk.X, padx=10, pady=(10, 5))

    Max_onS_option = tk.Radiobutton(master=right_frame, text='Max_onS', variable=boundary_var, value=2, indicator=0, background="light blue", command=select_boundary)
    Max_onS_option.pack(fill=tk.X, padx=10, pady=(10))


    # Radiobuttons for absolute value for the wave ('yes'/'no')
    absolute_var = tk.IntVar(master=right_frame, value=1)

    absolute_label = tk.Label(master=right_frame, text='Amplitude\n Absolute Value:')
    absolute_label.pack(fill=tk.X, padx=10, pady=(30, 0))

    no_absolute_option = tk.Radiobutton(master=right_frame, text='no', variable=absolute_var, value=1, indicator=0, background="light blue", command=select_absolute)
    no_absolute_option.pack(fill=tk.X, padx=10, pady=(10, 5))

    yes_absolute_option = tk.Radiobutton(master=right_frame, text='yes', variable=absolute_var, value=2, indicator=0, background="light blue", command=select_absolute)
    yes_absolute_option.pack(fill=tk.X, padx=10, pady=(10, 5))


    # Option to change the wave speed for wave animation
    speed_var = tk.StringVar(master=right_frame, value='0.75')

    speed_label = tk.Label(master=right_frame, text='Wave Speed:')
    speed_label.pack(fill=tk.X, padx=10, pady=(30, 0))

    speed_entry = tk.Entry(master=right_frame, textvariable=speed_var, width=10)
    speed_entry.pack(padx=10, pady=(10, 5))

    speed_button = tk.Button(master=right_frame, text='Apply Speed', command=change_speed, state='disabled')
    speed_button.pack(padx=10, pady=(10, 5))


    # Run the Tkinter main loop
    root.mainloop()




if __name__=='__main__':
    main()
