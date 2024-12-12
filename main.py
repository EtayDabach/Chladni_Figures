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

# Grid for rectangle membrane with boundary condition 0_onS
x0 = np.linspace(0,1,1000)
y0 = np.linspace(0,1,1000)
xv0, yv0 = np.meshgrid(x0,y0)

# Grid for rectangle membrane with boundary condition Max_onS
x = np.linspace(-0.5,0.5,1000)
y = np.linspace(-0.5,0.5,1000)
xv, yv = np.meshgrid(x,y) 

# Rectangle membrane dimensions
Lx = 1
Ly = 1

# Grid for circular membrane (symmetric and general) with boundary condition 0_onS
rj = np.linspace(0,1,1000) # Radius=1
thetaj= np.linspace(0, 2*np.pi, 1000)
rjv, thetajv = np.meshgrid(rj, thetaj)



# Initialize figure and axes based on membrane shape and animation type (plane or particles, but for now only particles)
class Membrane():
    def __init__(self, mtype:str, boundary='0_onS') -> None:
        self.type = mtype
        self.boundary = boundary
        # self.amplitude = 0
        self.shape()

    def shape(self) -> None:
        # Setting the initial parameters for
        if self.type == 'rect':
            self.fig, self.ax = plt.subplots(figsize=(8,6))
            if self.boundary == '0_onS':
                self.amplitude = wave_amp(x_vals=xv0, y_vals=yv0, n=0, m=0, L_x=Lx, L_y=Ly, boundary=self.boundary, mtype=self.type)
                self.ax.set_xlim(0, Lx)
                self.ax.set_ylim(0, Ly)
            elif self.boundary == 'Max_onS': # only odd modes for real resulst
                self.amplitude = wave_amp(x_vals=xv0, y_vals=yv0, n=0, m=0, L_x=Lx, L_y=Ly, boundary=self.boundary, mtype=self.type)
                self.ax.set_xlim(-Lx/2, Lx/2)
                self.ax.set_ylim(-Ly/2, Ly/2)

        elif 'circ' in self.type:
            self.fig, self.ax = plt.subplots(subplot_kw={'projection':'polar'}, figsize=(8,6))
            self.amplitude = wave_amp(x_vals=rjv, y_vals=thetajv, n=0, m=1, boundary='0_onS', mtype=self.type)
            self.ax.set_ylim(0, 1.05) # radius value + 0.05

        self.dot, = plt.plot([],[], 'o', ms=2, color='#cbbd93') # sand color #cbbd93
        self.ax.set_facecolor('dimgrey')    
        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)









if __name__=='__main__':
    pass
