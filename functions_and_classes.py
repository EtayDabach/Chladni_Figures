import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import scipy.special as sps
from matplotlib import animation
from mpl_toolkits import mplot3d
import ipywidgets as widg
from IPython.display import HTML
import IPython
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import tkinter as tk



def wave_amp(x_vals:np.ndarray, y_vals:np.ndarray, n:int, m:int, L_x=1, L_y=1, boundary='0_onS', mtype='rect' , abso='yes') -> np.ndarray:
    """
    Calculate the amplitude of the wave in a rectangular or circular membrane. 
    For rectangular membrane (mtype='rect') with boundary condition of maximum amplitude (boundary='Max_onS'), the modes n and m must be odd numbers to get the wanted result.
    For circular membrane with axisymmetric ,rotational symmetry, (mtype='circ_sym') n mode is always n=0 to get the bessel function of first kind in order 0 J_0.

    Args:
        x_vals (np.ndarray): Coordinates of x values for rectangular membrane or radius in polar coordinates for circular membrane.
        y_vals (np.ndarray): Coordinates of y values for rectangular membrane or theta in polar coordinates for circular membrane.
        n (int): The n mode of the wave in the x axis of the membrane (for circular membrane its reflected as the n order of the first kind bessel function J_n).
        m (int): The m mode of the wave in the y axis of the membrane (for circular membrane its reflected as the m zero of the J_n function, k_{n,m}).
        L_x (int, optional): Length of the x side for rectangular membrane (not affecting circular). Defaults to 1.
        L_y (int, optional): Length of the y side for rectangular membrane (not affecting circular). Defaults to 1.
        boundary (str, optional): Boundary condition of the membrane. For rectangular you can use '0_onS' or 'Max_onS', for circular only '0_onS'. Defaults to '0_onS'.
        mtype (str, optional): Define the type of the membrane as 'rect' , 'circ_sym' or 'circ_gen'. For 'circ_sym' the n mode value is always 0. Defaults to 'rect'.
        abso (str, optional): Option to get absolute value for the amplitude as 'yes' or 'no'. Defaults to 'yes'.

    Returns:
        np.ndarray: Returns an array with the values of the amplitudes.
    """
    if (mtype == 'rect') and (boundary == '0_onS'):
        amp = (np.sin(n*np.pi*x_vals / L_x) * np.sin(m*np.pi*y_vals / L_y))
    elif (mtype == 'rect') and (boundary == 'Max_onS'): 
        amp =  ((np.sin(n*np.pi*x_vals / L_x) * np.sin(m*np.pi*y_vals / L_y)) - (np.sin(m*np.pi*x_vals / L_x) * np.sin(n*np.pi*y_vals / L_y)))
    elif (mtype == 'circ_sym') and (boundary == '0_onS'): # symmertric circle (n=0)
        k_0m = sps.jn_zeros(0,m)[m-1]
        amp =  (sps.jv(0, k_0m*x_vals))
    elif (mtype == 'circ_gen') and (boundary == '0_onS'): # general circle (n,m)
        k_nm = sps.jn_zeros(n,m)[m-1]
        amp =  (sps.jv(n, k_nm*x_vals)*(np.cos(n*y_vals) + np.sin(n*y_vals)))
    
    if abso == 'no':
        return amp
    else:
        return abs(amp)



def time_evolution(omega:float, t:float) -> float:
    """
    Function for time evolution of the wave ,T(t), using the method of sepration of variables for both rectangular and circular membrane.
    The full function contains sum of sine and cosine, but for in this function its only contains the sine (initial condition in t=0 is amplitude=0).

        T(t)=a_{n,m}*cos(omega_{n,m}t) + b_{n,m}*sin(omega_{n,m}t)

    Args:
        omega (float): The angular velocity 
        t (float): Time

    Returns:
        float: Returns sine of omega*t
    """
    return np.sin(omega * t)



def rect_omega(n:int, m:int, L_x=1, L_y=1, speed_of_sound=1.0) -> float:
    """
    Calculate the angular velocity of the wave in rectangular membrane based on the modes n and m and the shape L_x and L_y.

    Args:
        n (int): The n mode of the wave in the x axis of the membrane
        m (int): The m mode of the wave in the y axis of the membrane
        L_x (int, optional): Length of the x side for rectangular membrane. Defaults to 1.
        L_y (int, optional): Length of the y side for rectangular membrane. Defaults to 1.
        speed_of_sound (float, optional): The speed of sound (wave speed) in the membrane. Defaults to 1.0.

    Returns:
        float: Returns the angular velocity, omega, in rectangular membrane.
    """
    return speed_of_sound * np.sqrt((n/L_x)**2 + (m/L_y)**2)



def circ_omega(n:int, m:int, radius=1.0, speed_of_sound=1.0) -> float:
    """
    Calculate the angular velocity of the wave in circular membrane based on the modes n and m and the radius.

    Args:
        n (int): The n mode reflected as the n order of the first kind bessel function J_n
        m (int): The m mode reflected as the m 'th zero of the J_n function, k_{n,m}
        radius (float, optional): Radius of the circular membrane. Defaults to 1.0.
        speed_of_sound (float, optional): The speed of sound (wave speed) in the membrane. Defaults to 1.0.

    Returns:
        float: Returns the angular velocity, omega, in circular membrane.
    """
    k_nm = sps.jn_zeros(n,m)[m-1]
    return (k_nm / radius) * speed_of_sound



class Particles:
    """
    Generate (points) particles from uniform distribution with random location on the membrane.
    Deafult membrane is 'rect' type.

    Attributes
    ----------
    amplitude : function
        A function that calculate the amplitude in the membrane.
    
    num_points : int
        Number of points (particles) to generate.

    a : int, b : int
        Limits to generates points for rectangular membrane. Needs to consider the size of the membrane (L_x and L_y) for better results.
        For now its only works for SQUARE membranes (L_x=L_y).
    
    ttype : str
        Type of the membrane as 'rect', 'circ_sym' or 'circ_gen'.
    
    delta : float
        Control the size of the step for the particles.
    
    prev_points : ndarray or int
        Previous position of the points

    """
    def __init__(self, amplitude=wave_amp, a=0, b=1, ttype='rect', num_points=10000, delta=0.05) -> None:
        """
        Constructs all the necessary attributes for the Particles object.

        Args:
            amplitude (function, optional): A function that calculate the amplitude in the membrane.. Defaults to wave_amp.
            a (int, optional): Minimum limit to generates points for rectangular membrane. Defaults to 0.
            b (int, optional): Maximum limits to generates points for rectangular membrane. Defaults to 1.
            ttype (str, optional): Type of the membrane as 'rect' or 'circ'. Defaults to 'rect'.
            num_points (int, optional): Number of points (particles) to generate. Defaults to 10000.
            delta (float, optional): Control the size of the step for the particles. Defaults to 0.05.
        """
        self.num_points = num_points
        self.amplitude = amplitude
        self.type=ttype
        self.prev_points = 0
        if 'circ' in self.type:
        # if self.type == 'circ':
            self.delta = delta/10
        else:
            self.delta = delta

        # Generate points based on the membrane shape
        if 'circ' in self.type:
        # if self.type == 'circ':
            polar_r_points = np.random.uniform(0, 1.0, self.num_points)
            polar_angle_points = np.random.uniform(0, 2*np.pi,self.num_points)
            self.points = np.array([polar_r_points,polar_angle_points])
        else:
            self.points = np.random.uniform(a, b, size=(2, self.num_points))

    
    def step(self, **amplitude_params) -> None:
        """
        Evolving all particle positions, dr, by generating a random direction angle from uniform distribution and multiply it by the normalized wave amplitude and constraining constant.
        Keeping track on previous locations to generate smooth animation.
        """
        n_mode = amplitude_params['n']
        m_mode = amplitude_params['m']

        # Normalization for the step based on the membrane shape
        if 'circ' in self.type:
        # if self.type == 'circ':
            k_m1n = sps.jn_zeros(n_mode, m_mode)[m_mode-1]
            norm = 1/2 * (sps.jv(n_mode+1, k_m1n))**2
        else:
            norm = 2

        angles = np.random.uniform(0, 2*np.pi, size=self.num_points)
        dr = self.delta * np.array([np.cos(angles), np.sin(angles)]) * self.amplitude(*self.points, **amplitude_params) / norm
        self.prev_points = np.copy(self.points)
        self.points += dr


