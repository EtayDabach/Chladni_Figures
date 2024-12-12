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
x0 = np.linspace(0, Lx, 1000)
y0 = np.linspace(0, Ly, 1000)
xv0, yv0 = np.meshgrid(x0, y0)

# Grid for rectangle membrane with boundary condition Max_onS
x = np.linspace(-Lx/2, Ly/2, 1000)
y = np.linspace(-Lx/2, Ly/2, 1000)
xv, yv = np.meshgrid(x, y) 


# Grid for circular membrane (symmetric and general) with boundary condition 0_onS
rj = np.linspace(0, radius, 1000) # radius=1
thetaj= np.linspace(0, 2*np.pi, 1000)
rjv, thetajv = np.meshgrid(rj, thetaj)



# Initialize figure and axes based on membrane shape and animation type (plane or particles, but for now only particles)
class Membrane():
    """_summary_
    """
    def __init__(self, mtype:str, boundary='0_onS') -> None:
        self.type = mtype
        self.boundary = boundary
        self.resolution = 1 # for slider intervals
        self.a = 0 # for uniform distribution
        self.b = 1

        # Setting initial amplitude in the membrane
        if self.type == 'rect':
            if self.boundary == '0_onS':
                self.amplitude = wave_amp(x_vals=xv0, y_vals=yv0, n=0, m=0, L_x=Lx, L_y=Ly, boundary=self.boundary, mtype=self.type)
                self.a = 0 # for uniform distribution
                self.b = Lx
            else:
                self.amplitude = wave_amp(x_vals=xv, y_vals=yv, n=0, m=0, L_x=Lx, L_y=Ly, boundary=self.boundary, mtype=self.type)
                self.resolution = 2 # for slider intervals to get only odd numbers
                self.a = -Lx/2 # for uniform distribution
                self.b = Lx/2

        elif 'circ' in self.type:
            self.amplitude = wave_amp(x_vals=rjv, y_vals=thetajv, n=0, m=1, boundary='0_onS', mtype=self.type)
            self.a = 0 # for uniform distribution
            self.b = 2*np.pi
        else:
            print("Please select one of the following options: 'rect', 'circ_sym' or 'circ_gen'.")
            self.amplitude = 0

        self.shape_for_particles()

    def shape_for_particles(self) -> None:
        # Setting the initial parameters for
        if self.type == 'rect':
            self.fig, self.ax = plt.subplots(figsize=(8,6))
            if self.boundary == '0_onS':
                self.ax.set_xlim(0, Lx)
                self.ax.set_ylim(0, Ly)
            elif self.boundary == 'Max_onS': # only odd modes for true resulst
                self.ax.set_xlim(-Lx/2, Lx/2)
                self.ax.set_ylim(-Ly/2, Ly/2)

        elif 'circ' in self.type:
            self.fig, self.ax = plt.subplots(subplot_kw={'projection':'polar'}, figsize=(8,6))
            self.ax.set_ylim(0, radius + 0.05) # radius value + 0.05 for better visualization

        self.dot, = plt.plot([],[], 'o', ms=2, color='#cbbd93') # sand color #cbbd93
        self.ax.set_facecolor('dimgrey')    
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)




def main():
    # Initialize the membrane
    membrane = Membrane(mtype='circ_sym', boundary='0_onS') # for mtype: 'rect', 'circ_sym' or 'circ_gen'.
    print(membrane.boundary)
    print(membrane.type)

    # Update function for the animation,
    def animate(i):
        global ensemble
        n_mode = slider_n.get()
        m_mode = slider_m.get()
        if i%5 == 0:
            ensemble.step(L_x=Lx, L_y=Ly, n=n_mode, m=m_mode, boundary=membrane.boundary, mtype=membrane.type, abso='yes')

        points = ensemble.prev_points + (i%5)/5 * (ensemble.points - ensemble.prev_points)
        membrane.dot.set_data(*points)
    
    # Create function for reset animation
    def create_particles():
        global ensemble
        ensemble = Particles(amplitude=wave_amp, a=membrane.a, b=membrane.b, ttype=membrane.type, num_points=10000, delta=0.1)
    
    # Initialization function for start animation
    def start_animation():
        global ani
        ani = animation.FuncAnimation(membrane.fig, animate, frames=200, interval=20, repeat=True)
        canvas.draw_idle()  # ensure the canvas updates

    # Stop function for stop animation
    def stop_animation():
        global ani
        if 'ani' in globals():
            ani.event_source.stop()  # stop the animation


    # Tkinter GUI setup
    root = tk.Tk()
    root.title("Interactive Chladni Figures")
    root.minsize(700, 700)

    create_particles() # initialize the particles
    print(ensemble.points)

    canvas = FigureCanvasTkAgg(membrane.fig, master=root)
    canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=1)

    # Control panel
    frame = tk.Frame(root)
    frame.pack(side=tk.BOTTOM)

    # Sliders for n and m
    slider_n = tk.Scale(frame, from_=1, to=15, resolution=membrane.resolution, label="n", orient="horizontal")
    slider_n.set(1)
    slider_n.pack(side=tk.LEFT, padx=10)

    slider_m = tk.Scale(frame, from_=1, to=15, resolution=membrane.resolution, label="m", orient="horizontal")
    slider_m.set(3)
    slider_m.pack(side=tk.LEFT, padx=10)

    # Buttons for start, stop and reset animation
    start_button = tk.Button(frame, text="Start Animation", command=start_animation)
    start_button.pack(side=tk.LEFT, padx=10)

    stop_button = tk.Button(frame, text="Stop Animation", command=stop_animation)
    stop_button.pack(side=tk.LEFT, padx=10)

    reset_button = tk.Button(frame, text="Reset Animation", command=create_particles) # , command=create_particles
    reset_button.pack(side=tk.LEFT, padx=10)

    # Run the Tkinter main loop
    root.mainloop()




if __name__=='__main__':
    main()
