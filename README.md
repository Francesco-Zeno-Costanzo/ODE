# ODE
These codes are examples of solving ordinary differential equations using the finite difference method, the method of runge kutta of order 4, scipy's "odeint" function and other methods.
Some codes also feature an animation of the solution.
The "neōn_katalogos" file (as the name might suggest) is a list containing several examples of implicit, explicit and symplectic ode-solving algorithms; all of them, for simplicity are applied to the harmonic oscillator. Some of them are also in the file ode.c

We briefly present the various methods used, each of these, unless otherwise specified, is contained as an example in the file neōn katalogos:

## Initiaal value problem
### Euler

The former are Euler's methods, classical, implicit and semi-implicit (the latter symplectic).
If we have the following:

$$
\begin{cases}
\displaystyle \frac{dy}{dt} = f(t,y) \\
y(t_0) = y_0
\end{cases}
$$

the explicit method says:

$$
y_{k+1} = y_k + \Delta t \hspace{1 mm} f(t_k, y_k)
$$

while the impicit:

$$
y_{k+1} = y_k + \Delta t \hspace{1 mm} f(t_{k+1}, y_{k+1})
$$

If we insted consider an hamiltonian system:

$$
\begin{cases}
\displaystyle \frac{dx}{dt} = f(t,v) \\
\displaystyle \frac{dv}{dt} = g(t,x)
\end{cases}
$$

We can use the semi-impicit method:

$$
\begin{aligned}
v_{k+1} &= v_k + \Delta t \hspace{1 mm} g(t_k, x_k) \\
x_{k+1} &= x_k + \Delta t \hspace{1 mm} f(t_k, v_{k+1})
\end{aligned}
$$

it is easy to verify the advantage of the latter method compared to the previous two as it is associated with a canonical transformation.

### Mid-point

Another first order method and symplectic integrator is the midpoint method that we present in its implicit and explicit formulation:
implicit:

$$
y_{k+1} = y_k + \Delta t \hspace{1 mm}
f \left(
t_k + \frac{\Delta t}{2},
\frac{y_k + y_{k+1}}{2}
\right)
$$

explicit:

$$
y_{k+1} = y_k + \Delta t \hspace{1 mm}
f \left(
t_k + \frac{\Delta t}{2},
y_k + \frac{\Delta t}{2} f(t_k, y_k)
\right)
$$

### Velocity verlet 

Another algorithm is verlet velocity integration, a symplectic method:

$$
\begin{aligned}
x_{k+1} &= x_k + v_k \Delta t + \frac{1}{2} a_k \Delta t^2 \\
v_{k+1} &= v_k + \frac{\Delta t}{2} \left( a_k + a_{k+1} \right)
\end{aligned}
$$

### Runge-Kutta 4

A method that could not miss the is :runge kutta of 4 th order:

$$ y_{k+1} = y_k + \frac{\Delta t}{6} \left( k_1 + 2k_2 + 2k_3 + k_4 \right)
$$

where:

$$
\begin{aligned}
k_1 &= f(t_k, y_k) \\
k_2 &= f \left(t_k + \frac{\Delta t}{2}, y_k + \frac{\Delta t}{2} k_1\right) \\
k_3 &= f \left(t_k + \frac{\Delta t}{2}, y_k + \frac{\Delta t}{2} k_2\right) \\
k_4 &= f(t_k + \Delta t, y_k + \Delta t \hspace{1 mm} k_3)
\end{aligned}
$$

if f=f(t) this method became Cavalieri-Simpson rule

### Yoshida4

If we wanted a symplectic integrator of higher order we can use Yoshida4:

$$
\begin{aligned}
v_{k+1} &= v_k + d_i \hspace{1 mm} a(x_k) \hspace{1 mm} \Delta t \\
x_{k+1} &= x_k + c_i \hspace{1 mm} v_{k+1} \hspace{1 mm} \Delta t
\end{aligned}
$$


the value of the coefficients d_i and c_i are in the code

### Prediction and correction

We can also use multiple integrators through the prediction and correction method.
In this code it is implemented using the classic euler and the trapezoidal rule.
With classical Euler we do the prevision and then correct them with trapezoidal rule:

$$
\begin{aligned}
\text{Predictor:}& \qquad
\overline{y}_{k+1} = y_k + \Delta t \hspace{1 mm} f(t_k, y_k)
\text{Corrector:}& \qquad
y_{k+1} = y_k + \frac{\Delta t}{2} \left[ f(t_k, y_k) + f(t_{k+1}, \overline{y}_{k+1})\right]
\end{aligned}
$$



### Adams-Bashforth-Moulton, order 4

The following is a fourth-order predictor-corrector method using an explicit Adams-Bashforth scheme as a predictor and an implicit Adams-Moulton scheme as a corrector. To get started, three points of the solution are needed, which are calculated using an RK4

$$
\begin{aligned}
\text{Predictor:}& \quad
\overline{x}_{i+1} = x_i + \frac{h}{24}\left( 55 f_i - 59 f_{i-1} + 37 f_{i-2} - 9 f_{i-3}\right)
\text{Corrector:}& \quad
x_{i+1} = x_i + \frac{h}{24} \left( 9 f_{i+1} + 19 f_i - 5 f_{i-1} + f_{i-2} \right)
\end{aligned}
$$


## Boundary problem
### Shooting

A technique that is not used in the neōn katalogos code is the shooting method, because it was made for  boundary value problems and not for initial value problems.
It allows you to solve a boundary value problem by reducing it to the system of an initial value problem as we can see in shooting.py.
In this code we use the first derivative in t = 0  as a parameter that is not known and we try to find the value for which the solution takes the correct value on the boundary.
(we remember that in general the solution isn't unique for boundary problems)

If we have the following boundary value problem:

$$
\begin{cases}
\ddot{y}(t) = f(t, y(t), \dot{y}(t)) \\
y(t_0) = y_0 \\
y(t_1) = y_1
\end{cases}
$$

We can rewrite it as an initial value problem:

$$
\begin{cases}
\ddot{y}(t) = f(t, y(t), \dot{y}(t)) \\
y(t_0) = y_0 \\
\dot{y}(t_0) = s
\end{cases}
$$

So now the problem become to find the solution of the following equation:

$$
F(s) = y(t_1; s) - y_1 = 0
$$

In some particular case, as in Schrodinger's equation where the initial condition are known, if the potential is even, we can use the energy as parameter.
Two examples are in buca quadrata e osc.arm-s (see repository: quantum-mechanics)

### Relaxation





## Adaptive integrator

The adaptive.py code is an example of an integrator with an adaptive step (Runge–Kutta–Fehlberg method and Cash-Karp Runga-Kutta Method), i.e. the integration step changes during the simulation. 

### Runge–Kutta–Fehlberg
If we have an initial value problem:

$$
\begin{cases}
\displaystyle \frac{dy}{dt} = f(t,y) \\
y(t_0) = y_0
\end{cases}
$$

then we write the solution as:

$$
y_{k+1} = y_k + \sum_{i=1}^{6} CH_i k_i
$$

where:

$$
\begin{aligned}
k_1 &= h f(t + A_0 h, y) \\
k_2 &= h f(t + A_1 h, y + B_{11} k_1) \\
k_3 &= h f(t + A_2 h, y + B_{21} k_1 + B_{22} k_2) \\
k_4 &= h f(t + A_3 h, y + B_{31} k_1 + B_{32} k_2 + B_{33} k_3) \\
k_5 &= h f(t + A_4 h, y + B_{41} k_1 + B_{42} k_2 + B_{43} k_3 + B_{44} k_4) \\
k_6 &= h f(t + A_5 h, y + B_{51} k_1 + B_{52} k_2 + B_{53} k_3 + B_{54} k_4 + B_{55} k_5)
\end{aligned}
$$

The extimention of the truncation error is:

$$
TE = \left| \sum_{i=1}^{6} CT_i k_i \right|
$$

and the update of the update of the step is done in the following way:

$$
h_{\text{new}} = 0.9 h \left( \frac{\varepsilon}{TE} \right)^{1/5}
$$

The coefficients are reported in the code.

### 
Cash-Karp Runga-Kutta