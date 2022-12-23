# ODE
These codes are examples of solving ordinary differential equations using the finite difference method, the method of runge kutta of order 4, scipy's "odeint" function and other methods.
Some codes also feature an animation of the solution.
The "neōn_katalogos" file (as the name might suggest) is a list containing several examples of implicit, explicit and symplectic ode-solving algorithms; all of them, for simplicity are applied to the harmonic oscillator. Some of them are also in the file ode.c

We briefly present the various methods used, each of these, unless otherwise specified, is contained as an example in the file neōn katalogos:

## Euler

The former are Euler's methods, classical, implicit and semi-implicit (the latter symplectic)


<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;\begin{cases}&space;\frac{dy}{dt}=f(t,&space;y)\\&space;y(t_0)=y_0&space;\end{cases}&space;\\&space;\text{classical&space;euler&space;method:}\\&space;y_{k&plus;1}=y_k&plus;dtf(t_k,&space;y_k)&space;\\&space;\\&space;\text{implicit&space;euler&space;method:}\\&space;y_{k&plus;1}=y_k&plus;dtf(t_{k&plus;1},&space;y_{k&plus;1})\\&space;\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;\begin{cases}&space;\frac{dy}{dt}=f(t,&space;y)\\&space;y(t_0)=y_0&space;\end{cases}&space;\\&space;\text{classical&space;euler&space;method:}\\&space;y_{k&plus;1}=y_k&plus;dtf(t_k,&space;y_k)&space;\\&space;\\&space;\text{implicit&space;euler&space;method:}\\&space;y_{k&plus;1}=y_k&plus;dtf(t_{k&plus;1},&space;y_{k&plus;1})\\&space;\\" title="\\ \begin{cases} \frac{dy}{dt}=f(t, y)\\ y(t_0)=y_0 \end{cases} \\ \text{classical euler method:}\\ y_{k+1}=y_k+dtf(t_k, y_k) \\ \\ \text{implicit euler method:}\\ y_{k+1}=y_k+dtf(t_{k+1}, y_{k+1})\\ \\" /></a>



<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;\text{the&space;semi-implicit&space;euler&space;method&space;is&space;good}\\&space;\text{when&space;dealing&space;with&space;the&space;Hamilton's&space;equations}\\&space;\text{which&space;are&space;of&space;the&space;form:}\\&space;\\&space;\begin{cases}&space;\frac{dx}{dt}=f(t,&space;v)\\&space;\frac{dv}{dt}=g(t,&space;x)\\&space;\end{cases}&space;\\&space;\text{so&space;we&space;have:}\\&space;v_{k&plus;1}=v_k&plus;dtg(t_{k},&space;x_{k})\\&space;x_{k&plus;1}=x_k&plus;dtf(t_{k},&space;v_{k&plus;1})\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;\text{the&space;semi-implicit&space;euler&space;method&space;is&space;good}\\&space;\text{when&space;dealing&space;with&space;the&space;Hamilton's&space;equations}\\&space;\text{which&space;are&space;of&space;the&space;form:}\\&space;\\&space;\begin{cases}&space;\frac{dx}{dt}=f(t,&space;v)\\&space;\frac{dv}{dt}=g(t,&space;x)\\&space;\end{cases}&space;\\&space;\text{so&space;we&space;have:}\\&space;v_{k&plus;1}=v_k&plus;dtg(t_{k},&space;x_{k})\\&space;x_{k&plus;1}=x_k&plus;dtf(t_{k},&space;v_{k&plus;1})\\" title="\\ \text{the semi-implicit euler method is good}\\ \text{when dealing with the Hamilton's equations}\\ \text{which are of the form:}\\ \\ \begin{cases} \frac{dx}{dt}=f(t, v)\\ \frac{dv}{dt}=g(t, x)\\ \end{cases} \\ \text{so we have:}\\ v_{k+1}=v_k+dtg(t_{k}, x_{k})\\ x_{k+1}=x_k+dtf(t_{k}, v_{k+1})\\" /></a>


it is easy to verify the advantage of the latter method compared to the previous two as it is associated with a canonical transformation.

## Mid-point

Another first order method and symplectic integrator is the midpoint method that we present in its implicit and explicit formulation:
implicit:

<a href="https://www.codecogs.com/eqnedit.php?latex=y_{k&plus;1}=y_k&plus;dtf\Biggl(t_{k}&plus;\frac{dt}{2},\frac{1}{2}(y_k&space;&plus;&space;y_{k&plus;1})\Biggr)\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{k&plus;1}=y_k&plus;dtf\Biggl(t_{k}&plus;\frac{dt}{2},\frac{1}{2}(y_k&space;&plus;&space;y_{k&plus;1})\Biggr)\\" title="y_{k+1}=y_k+dtf\Biggl(t_{k}+\frac{dt}{2},\frac{1}{2}(y_k + y_{k+1})\Biggr)\\" /></a>

explicit:

<a href="https://www.codecogs.com/eqnedit.php?latex=y_{k&plus;1}=y_k&plus;dtf\Biggl(t_{k}&plus;\frac{dt}{2},&space;y_k&space;&plus;&space;\frac{dt}{2}f(t_k,&space;y_k)\Biggr)\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?y_{k&plus;1}=y_k&plus;dtf\Biggl(t_{k}&plus;\frac{dt}{2},&space;y_k&space;&plus;&space;\frac{dt}{2}f(t_k,&space;y_k)\Biggr)\\" title="y_{k+1}=y_k+dtf\Biggl(t_{k}+\frac{dt}{2}, y_k + \frac{dt}{2}f(t_k, y_k)\Biggr)\\" /></a>

## velocity verlet 

Another algorithm is verlet velocity integration, a symplectic method:

<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;x_{k&plus;1}=x_{k}&plus;v_{k}dt&space;&plus;a_{k}&space;\frac{dt^2}{2}\\&space;v_{k&plus;1}=v_{k}&plus;\frac{dt}{2}(a_{k}&plus;a_{k&plus;1})" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;x_{k&plus;1}=x_{k}&plus;v_{k}dt&space;&plus;a_{k}&space;\frac{dt^2}{2}\\&space;v_{k&plus;1}=v_{k}&plus;\frac{dt}{2}(a_{k}&plus;a_{k&plus;1})" title="\\ x_{k+1}=x_{k}+v_{k}dt +a_{k} \frac{dt^2}{2}\\ v_{k+1}=v_{k}+\frac{dt}{2}(a_{k}+a_{k+1})" /></a>

## Runge-Kutta 4

A method that could not miss the is :runge kutta of 4 th order:


<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;\begin{cases}&space;\dot{y}(t)=f(t,y(t))\\&space;y(t_0)=y_0&space;\end{cases}\\&space;y_{k&plus;1}&space;=&space;y_{k}&space;&plus;&space;\tfrac{h}{6}\left(k_1&space;&plus;&space;2k_2&space;&plus;&space;2k_3&space;&plus;&space;k_4&space;\right)\\&space;t_{k&plus;1}&space;=&space;t_k&space;&plus;&space;dt&space;\\&space;\text{where}&space;\hspace{2&space;mm}&space;k_1,&space;k_2,&space;k_3,&space;k_4&space;\hspace{2&space;mm}&space;\text{are:}&space;\\&space;k_1&space;=&space;f(t_k,&space;y_n)&space;\\&space;k_2&space;=&space;f(t_k&space;&plus;&space;\tfrac{dt}{2},&space;y_n&space;&plus;&space;\tfrac{1}{2}&space;k_1&space;dt)&space;\\&space;k_3&space;=&space;f(t_k&space;&plus;&space;\tfrac{dt}{2},&space;y_n&space;&plus;&space;\tfrac{1}{2}&space;k_2&space;dt)&space;\\&space;k_4&space;=&space;f(t_k&space;&plus;&space;dt,&space;y_n&space;&plus;&space;k_3&space;dt)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;\begin{cases}&space;\dot{y}(t)=f(t,y(t))\\&space;y(t_0)=y_0&space;\end{cases}\\&space;y_{k&plus;1}&space;=&space;y_{k}&space;&plus;&space;\tfrac{h}{6}\left(k_1&space;&plus;&space;2k_2&space;&plus;&space;2k_3&space;&plus;&space;k_4&space;\right)\\&space;t_{k&plus;1}&space;=&space;t_k&space;&plus;&space;dt&space;\\&space;\text{where}&space;\hspace{2&space;mm}&space;k_1,&space;k_2,&space;k_3,&space;k_4&space;\hspace{2&space;mm}&space;\text{are:}&space;\\&space;k_1&space;=&space;f(t_k,&space;y_n)&space;\\&space;k_2&space;=&space;f(t_k&space;&plus;&space;\tfrac{dt}{2},&space;y_n&space;&plus;&space;\tfrac{1}{2}&space;k_1&space;dt)&space;\\&space;k_3&space;=&space;f(t_k&space;&plus;&space;\tfrac{dt}{2},&space;y_n&space;&plus;&space;\tfrac{1}{2}&space;k_2&space;dt)&space;\\&space;k_4&space;=&space;f(t_k&space;&plus;&space;dt,&space;y_n&space;&plus;&space;k_3&space;dt)" title="\\ \begin{cases} \dot{y}(t)=f(t,y(t))\\ y(t_0)=y_0 \end{cases}\\ y_{k+1} = y_{k} + \tfrac{h}{6}\left(k_1 + 2k_2 + 2k_3 + k_4 \right)\\ t_{k+1} = t_k + dt \\ \text{where} \hspace{2 mm} k_1, k_2, k_3, k_4 \hspace{2 mm} \text{are:} \\ k_1 = f(t_k, y_n) \\ k_2 = f(t_k + \tfrac{dt}{2}, y_n + \tfrac{1}{2} k_1 dt) \\ k_3 = f(t_k + \tfrac{dt}{2}, y_n + \tfrac{1}{2} k_2 dt) \\ k_4 = f(t_k + dt, y_n + k_3 dt)" /></a>

if f=f(t) this method became Cavalieri-Simpson rule

## Yoshida4

If we wanted a symplectic integrator of higher order we can use Yoshida4:

<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;v_{k&plus;1}=v_{k}&plus;d_{i}a(x_{k})dt\\&space;x_{k&plus;1}=x_{k}&plus;c_{i}v_{k&plus;1}dt\\" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;v_{k&plus;1}=v_{k}&plus;d_{i}a(x_{k})dt\\&space;x_{k&plus;1}=x_{k}&plus;c_{i}v_{k&plus;1}dt\\" title="\\ v_{k+1}=v_{k}+d_{i}a(x_{k})dt\\ x_{k+1}=x_{k}+c_{i}v_{k+1}dt\\" /></a>


the value of the coefficients d_i and c_i are in the code

## prediction and correction

We can also use multiple integrators through  the prediction and correction method.
In this code it is implemented using the classic euler and the trapezoidal rule.
With classical Euler we do the prevision and then correct them with trapezoidal rule:


<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;\text{Predictor:}\\&space;\overline{y}_{k&plus;1}=y_k&plus;dtf(t_k,&space;y_k)\\&space;\text{Corrector:}\\&space;y_{k&plus;1}&space;=&space;y_k&space;&plus;&space;\frac{dt}{2}&space;\bigl(&space;f(t_k,&space;y_k)&space;&plus;&space;f(t_{k&plus;1},&space;\overline{y}_{k&plus;1})&space;\bigr)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;\text{Predictor:}\\&space;\overline{y}_{k&plus;1}=y_k&plus;dtf(t_k,&space;y_k)\\&space;\text{Corrector:}\\&space;y_{k&plus;1}&space;=&space;y_k&space;&plus;&space;\frac{dt}{2}&space;\bigl(&space;f(t_k,&space;y_k)&space;&plus;&space;f(t_{k&plus;1},&space;\overline{y}_{k&plus;1})&space;\bigr)" title="\\ \text{Predictor:}\\ \overline{y}_{k+1}=y_k+dtf(t_k, y_k)\\ \text{Corrector:}\\ y_{k+1} = y_k + \frac{dt}{2} \bigl( f(t_k, y_k) + f(t_{k+1}, \overline{y}_{k+1}) \bigr)" /></a>

## adams-bashforth-moulton, order 4

to do

## Boundary problem

A technique that is not used in the neōn katalogos code is the shooting method, because it was made for  boundary value problems and not for initial value problems.
It allows you to solve a boundary value problem by reducing it to the system of an initial value problem as we can see in shooting.py.
In this code we use the first derivative in t = 0  as a parameter that is not known and we try to find the value for which the solution takes the correct value on the boundary.
(we remember that in general the solution isn't unique for boundary problems)


<a href="https://www.codecogs.com/eqnedit.php?latex=\\&space;\begin{cases}&space;\ddot{y}(t)=f(t,&space;y(t),&space;\dot{y}(t))\\&space;y(t_0)=y_0\\&space;y(t_1)=y_1&space;\end{cases}\\&space;\begin{cases}&space;\ddot{y}(t)=f(t,&space;y(t),&space;\dot{y}(t))\\&space;y(t_0)=y_0\\&space;\dot{y}(t_0)=s&space;\end{cases}\\&space;\text{so&space;after&space;solving&space;the&space;second&space;system&space;with&space;different}&space;\hspace{1&space;mm}&space;s\\&space;\text{we&space;look&space;for}&space;\hspace{1&space;mm}&space;s&space;\hspace{1&space;mm}&space;\text{such&space;that&space;the&space;solution&space;with&space;that}\\&space;\hspace{1&space;mm}&space;s&space;\hspace{1&space;mm}&space;\text{&space;is&space;worth&space;in}\hspace{1&space;mm}&space;t_1&space;\hspace{1&space;mm}&space;\text{&space;as&space;much&space;as&space;we&space;want}\\&space;F(s)=y(t_1;s)-y_1=0" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\\&space;\begin{cases}&space;\ddot{y}(t)=f(t,&space;y(t),&space;\dot{y}(t))\\&space;y(t_0)=y_0\\&space;y(t_1)=y_1&space;\end{cases}\\&space;\begin{cases}&space;\ddot{y}(t)=f(t,&space;y(t),&space;\dot{y}(t))\\&space;y(t_0)=y_0\\&space;\dot{y}(t_0)=s&space;\end{cases}\\&space;\text{so&space;after&space;solving&space;the&space;second&space;system&space;with&space;different}&space;\hspace{1&space;mm}&space;s\\&space;\text{we&space;look&space;for}&space;\hspace{1&space;mm}&space;s&space;\hspace{1&space;mm}&space;\text{such&space;that&space;the&space;solution&space;with&space;that}\\&space;\hspace{1&space;mm}&space;s&space;\hspace{1&space;mm}&space;\text{&space;is&space;worth&space;in}\hspace{1&space;mm}&space;t_1&space;\hspace{1&space;mm}&space;\text{&space;as&space;much&space;as&space;we&space;want}\\&space;F(s)=y(t_1;s)-y_1=0" title="\\ \begin{cases} \ddot{y}(t)=f(t, y(t), \dot{y}(t))\\ y(t_0)=y_0\\ y(t_1)=y_1 \end{cases}\\ \begin{cases} \ddot{y}(t)=f(t, y(t), \dot{y}(t))\\ y(t_0)=y_0\\ \dot{y}(t_0)=s \end{cases}\\ \text{so after solving the second system with different} \hspace{1 mm} s\\ \text{we look for} \hspace{1 mm} s \hspace{1 mm} \text{such that the solution with that}\\ \hspace{1 mm} s \hspace{1 mm} \text{ is worth in}\hspace{1 mm} t_1 \hspace{1 mm} \text{ as much as we want}\\ F(s)=y(t_1;s)-y_1=0" /></a>


In some particular case, as in Schrodinger's equation where the initial condition are known, if the potential is even, we can use the energy as parameter.
Two examples are in buca quadrata e osc.arm-s

## Adaptive integrator

The RK4(5).py code is an example of an integrator with an adaptive step (Runge–Kutta–Fehlberg method), i.e. the integration step changes during the simulation. To test it, as always, the harmonic oscillator is used.
Brief explanation of the method:


<img src="https://latex.codecogs.com/svg.image?\\\begin{cases}\frac{dy}{dt}&space;=&space;f(t,&space;y(t))\\y(t_0)&space;=&space;y_0\\\end{cases}\\\begin{align*}k_1&=h&space;f(t&plus;A(0)h,&space;y)\\k_2&=h&space;f(t&plus;A(1)h,&space;y&plus;B1(1)k_1)\\k_3&=h&space;f(t&plus;A(2)h,&space;y&plus;B1(2)k_1&plus;B2(2)k_2&space;)\\k_4&=h&space;f(t&plus;A(3)h,&space;y&plus;B1(3)k_1&plus;B2(3)k_2&plus;B3(3)k_3&space;)\\k_5&=h&space;f(t&plus;A(4)h,&space;y&plus;B1(4)k_1&plus;B2(4)k_2&plus;B3(4)k_3&plus;B4(4)k_4&space;)\\k_6&=h&space;f(t&plus;A(5)h,&space;y&plus;B1(5)k_1&plus;B2(5)k_2&plus;B3(5)k_3&plus;B4(5)k_4&plus;B5(5)k_5)\end{align*}\\y_{k&plus;1}=y_k&plus;CH(0)k_1&plus;CH(1)k_2&plus;CH(2)k_3&plus;CH(3)k_4&plus;CH(4)k_5&plus;CH(5)k_6&space;\\\\\text{The&space;estimate&space;of&space;the&space;truncation&space;error&space;is:}\\TE=|CT(0)&space;k_1&space;&plus;CT(1)&space;k_2&plus;CT(2)k_3&plus;CT(3)&space;k_4&plus;CT(4)&space;k_5&plus;CT(5)&space;k_6|\\\\\text{At&space;the&space;completion&space;of&space;the&space;step,&space;a&space;new&space;stepsize&space;is&space;calculated:}\\h_{new}=0.9&space;&space;h&space;\left&space;(&space;\frac{\epsilon}{TE}&space;\right&space;)^{1/5}\\\text{If&space;TE&space;}&space;>&space;\epsilon,&space;\text{then&space;replace&space;h&space;with&space;}&space;h_{new}&space;\text{&space;and&space;repeat&space;the&space;step.&space;If&space;not&space;then&space;the&space;step&space;is&space;completed.}" title="\\\begin{cases}\frac{dy}{dt} = f(t, y(t))\\y(t_0) = y_0\\\end{cases}\\\begin{align*}k_1&=h f(t+A(0)h, y)\\k_2&=h f(t+A(1)h, y+B1(1)k_1)\\k_3&=h f(t+A(2)h, y+B1(2)k_1+B2(2)k_2 )\\k_4&=h f(t+A(3)h, y+B1(3)k_1+B2(3)k_2+B3(3)k_3 )\\k_5&=h f(t+A(4)h, y+B1(4)k_1+B2(4)k_2+B3(4)k_3+B4(4)k_4 )\\k_6&=h f(t+A(5)h, y+B1(5)k_1+B2(5)k_2+B3(5)k_3+B4(5)k_4+B5(5)k_5)\end{align*}\\y_{k+1}=y_k+CH(0)k_1+CH(1)k_2+CH(2)k_3+CH(3)k_4+CH(4)k_5+CH(5)k_6 \\\\\text{The estimate of the truncation error is:}\\TE=|CT(0) k_1 +CT(1) k_2+CT(2)k_3+CT(3) k_4+CT(4) k_5+CT(5) k_6|\\\\\text{At the completion of the step, a new stepsize is calculated:}\\h_{new}=0.9 h \left ( \frac{\epsilon}{TE} \right )^{1/5}\\\text{If TE } > \epsilon, \text{then replace h with } h_{new} \text{ and repeat the step. If not then the step is completed.}" />

The coefficients are reported in the code. It doesn't always work well unfortunately
