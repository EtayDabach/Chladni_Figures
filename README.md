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

For $u(x,y,t)=X(x)Y(y)T(t)$ and for boundary conditions of $u(x,y,t) = 0 ,\space \forall x,y \in \Sigma$ (fixed along all four edges) and inital conditions $u(x,y,t=0)=f(x,y)$ and $\partial_t{u(t=0)}=g(x,y)$, the amplitude is: $$\psi_{n,m} (x,y)=X(x)Y(y)=\sin{(\frac{n\pi}{L_x}x)}\sin{(\frac{m\pi}{L_y}y)}$$
and the frequencies are: 

$$f_{nm} = \frac{v}{2\pi}\sqrt{(\frac{n}{L_x})^2 + (\frac{m}{L_y})^2}, \space n,m = 1,2,3,4...$$

For different boundary conditions we get different amplitude. For example if we take the condition of maximum amplitude on $\Sigma$ we get: 

$$\psi_{n,m} (x,y)=\sin{(\frac{n\pi}{L_x}x)}\sin{(\frac{m\pi}{L_y}y)} - \sin{(\frac{m\pi}{L_x}x)}\sin{(\frac{n\pi}{L_y}y)}, \space n,m = 1,3,5,7...$$


<h2> Circular Membrane: </h2>

<h3> Rotational symmetry: </h3>

If we have rotational symmetry we get $u(r,\theta,t)=u(r,t)=R(r)T(t)$ in cylidrical coordinates, and for boundary conditions of $u(r=R,t) = 0,\space 0\leq r \leq R$ and inital conditions $u(t=0)=f(r)$ and $\partial_t{u(t=0)}=g(r)$ we get: 

$$R(r)=AJ_{0}(\frac{k_{0,m}}{R}r); \space R(r=R)=AJ_{0}(k_{0,m})=0 \space \Rightarrow k_{0,m} \space -roots  \space of J_0 $$

$J_0$ is the bessel function of first kind in order 0 and we get that $k_{0,m}$ are the roots of $J_0, \space m=1,2,3...$ $\\$
$\\$

<h3> General case: </h3>

$u(r,\theta,t)=R(r)\Theta(\theta)T(t)$ boundary conditions of $u(r=R,\theta,t) = 0,\space 0\leq r \leq R$ and inital conditions $u(t=0)=f(r,\theta)$ and $\partial_t{u(t=0)}=g(r,\theta)$:

$$\psi_{r,\theta}=R(r)\Theta(\theta)=AJ_{n}(\frac{k_{n,m}}{R}r)(a_{n}\cos{n\theta} + b_{n}\sin{n\theta}), \space m=1,2,3..; \space n=0,1,2,3....$$

<br>

The frequencies for both cases are: 

$$f_{n,m} = \frac{k_{n,m}}{2\pi R}*v$$

<hr>

<h2> GUI explanation (will be added soon with a short video/gif) </h2>
