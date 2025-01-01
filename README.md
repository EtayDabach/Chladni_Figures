<h1 style="text-align: center;"> Chladni Figures Simulation </h1>

2D wave equation: 

$$\frac{\partial^2 u}{\partial^2 t} = v^2 \nabla^2{u}$$

Solution is given by using method sepration of variable $u(x,y,t)=X(x)Y(y)T(t)$ or for polar coordiantes $u(r,\theta,t)=R(r)\Theta(\theta)T(t)$.

For both. the equation for the time , $T(t)$, is the same:

$$ T(t) = a_{n,m}\cos{(\omega_{n,m}t)} + b_{n,m}\sin{(\omega_{n,m}t)} $$

Given an initial condition at $t=t_{0}$ this equation can be solved. For $t_{0}=0$ we get:

$$ T(t) =  b_{n,m}\sin{(\omega_{n,m}t)} $$

I will choose $b_{n,m}=1$ and apply the normalization constant on the remaining functions (e.g $X(x)Y(y)$...), and eventually get:

$$ T(t) =  \sin{(\omega_{n,m}t)} $$

This will be applied in the waves animation along with the corresponding functions for the amplitude (and angles for circular membrane).

<h2> Rectangular Membrane: </h2>

For $u(x,y,t)=X(x)Y(y)T(t)$ and for boundary conditions of $u(x,y,t) = 0 ,\space \forall x,y \in \Sigma$ (fixed along all four edges) and inital conditions $u(x,y,t=0)=f(x,y)$ and $\partial_t{u(t=0)}=g(x,y)$, the amplitude is: $$\psi_{n,m} (x,y)=X(x)Y(y)=\sin{\left(\frac{n\pi}{L_x}x\right)} \sin{\left(\frac{m\pi}{L_y}y\right)}$$
and the frequencies are: 

$$f_{nm} = \frac{v}{2\pi}\sqrt{\left(\frac{n}{L_x}\right)^2 + \left(\frac{m}{L_y}\right)^2}, \space n,m = 1,2,3,4...$$

For different boundary conditions we get different amplitude. For example if we take the condition of maximum amplitude on $\Sigma$ we get: 

$$\psi_{n,m} (x,y)=\sin{\left(\frac{n\pi}{L_x}x\right)}\sin{\left(\frac{m\pi}{L_y}y\right)} - \sin{\left(\frac{m\pi}{L_x}x\right)}\sin{\left(\frac{n\pi}{L_y}y\right)}, \space n,m = 1,3,5,7...$$

Notice that $n \space \text{and} \space m$ are odd numbers only and for $n=m$ the amplitude is always zero.


<h2> Circular Membrane: </h2>

<h3> Rotational symmetry: </h3>

If we have rotational symmetry we get $u(r,\theta,t)=u(r,t)=R(r)T(t)$ in cylidrical coordinates, and for boundary conditions of $u(r=R,t) = 0,\space 0\leq r \leq R$ and inital conditions $u(t=0)=f(r)$ and $\partial_t{u(t=0)}=g(r)$ we get: 

$$R(r)=AJ_{0}\left(\frac{k_{0,m}}{R}r\right); \space R(r=R)=AJ_{0}(k_{0,m})=0 \space \Rightarrow k_{0,m} \space -\text{roots of} \space J_0 $$

$J_0$ is the bessel function of first kind in order 0 and we get that $k_{0,m}$ are the roots of $J_0, \space m=1,2,3...$

<h3> General case: </h3>

$u(r,\theta,t)=R(r)\Theta(\theta)T(t)$ boundary conditions of $u(r=R,\theta,t) = 0,\space 0\leq r \leq R$ and inital conditions $u(t=0)=f(r,\theta)$ and $\partial_t{u(t=0)}=g(r,\theta)$:

$$\psi_{r,\theta}=R(r)\Theta(\theta)=AJ_{n}\left(\frac{k_{n,m}}{R}r\right)(a_{n}\cos{n\theta} + b_{n}\sin{n\theta}), \space m=1,2,3..; \space n=0,1,2,3....$$

<br>

The frequencies for both cases are: 

$$f_{n,m} = \frac{k_{n,m}}{2\pi R}*v$$

<hr>


<h2> GUI explanation </h2>

![Chladni Figures GUI img](https://github.com/user-attachments/assets/a1c1575f-9dc1-4ad9-a7ba-527a2f5c50ba)

On the right side you can see all of the options for the simulation which allow you choose between two types of animation: Paticles and Waves for each of the membrane we discussed earlier, and on the buttom you can initiate the animation and change the values of the modes $n \space \text{and} \space m$.


<h3> Animation Options </h3>

Choose between Particles or Waves animation.

https://github.com/user-attachments/assets/660212b5-4994-4d19-963c-1035fec5c8b9


<h3> Shape Options </h3>

Choose between two membranes: Rectangular or Circular shape.

https://github.com/user-attachments/assets/815d1a5e-1d42-4e0c-b7ee-dde300895e9e

For circular shape you can choose between symmetrical or general case.


<h3> Boundary Options </h3>

Choose boundary condition for rectangular shape only (circular shape is always 0_onS).

https://github.com/user-attachments/assets/948cd69b-bdb1-44fd-b449-63d792b358d6


<h3> Particles Step (For Particles Animation Only) </h3>

Controls the particles step value by multiply it with the given input, making it easier to reveal the patterns in some occasions.

https://github.com/user-attachments/assets/7b012f09-606b-499f-bce5-e3601b99ec9b


<h3> Wave Speed (For Waves Animation Only) </h3>

Controls the wave speed in the membrane.

https://github.com/user-attachments/assets/5a594797-9de6-41d4-a9cd-1aa975137b4d


<h3> Amplitude Absolute Value </h3>


Allow you to get the absolute value of the amplitude. This option is mostly to reveal the patterns of which the particles allign from the particles animation in waves animation0.


https://github.com/user-attachments/assets/407e0555-ad74-46cd-b986-7e1da37da1a1






