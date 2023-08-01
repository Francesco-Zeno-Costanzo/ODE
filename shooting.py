'''
Boundary value problem
x''(t) = 3/2 x(t)
x(0) = x0
x(1) = x1
'''
import numpy as np
import  matplotlib.pyplot  as  plt

#============================================================================
# Itegration: Adams-Bashforth-Moulton predictor and corretor of order 4
#============================================================================

def AMB4(num_steps, t0, tf, f, init, args=()):
    """
    Integrator with Adams-Bashforth-Moulton
    predictor and corretor of order 4

    Parameters
    ----------
    num_steps : int
        number of point of solution
    t0 : float
        lower bound of integration
    tf : float
        upper bound of integration
    f : callable
        function to integrate, must accept vectorial input
    init : 1darray
        array of initial condition
    args : tuple, optional
        extra arguments to pass to f

    Return
    ------
    X : array, shape (num_steps + 1, len(init))
        solution of equation
    t : 1darray
        time
    """
    #time steps
    dt = tf/num_steps

    X = np.zeros((num_steps + 1, len(init))) #matrice delle soluzioni
    t = np.zeros(num_steps + 1)              #array dei tempi

    X[0, :] = init                           #condizioni iniziali
    t[0]    = t0

    #primi passi con runge kutta
    for i in range(3):
        xk1 = f(t[i], X[i, :], *args)
        xk2 = f(t[i] + dt/2, X[i, :] + xk1*dt/2, *args)
        xk3 = f(t[i] + dt/2, X[i, :] + xk2*dt/2, *args)
        xk4 = f(t[i] + dt, X[i, :] + xk3*dt, *args)
        X[i + 1, :] = X[i, :] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
        t[i + 1] = t[i] + dt

    # Adams-Bashforth-Moulton
    i = 3
    AB0 = f(t[i  ], X[i,   :], *args)
    AB1 = f(t[i-1], X[i-1, :], *args)
    AB2 = f(t[i-2], X[i-2, :], *args)
    AB3 = f(t[i-3], X[i-3, :], *args)

    for i in range(3,num_steps):
        #predico
        X[i + 1, :] = X[i, :] + dt/24*(55*AB0 - 59*AB1 + 37*AB2 - 9*AB3)
        t[i + 1] = t[i] + dt
        #correggo
        AB3 = AB2
        AB2 = AB1
        AB1 = AB0
        AB0 = f(t[i+1], X[i + 1, :], *args)

        X[i + 1, :] = X[i, :] + dt/24*(9*AB0 + 19*AB1 - 5*AB2 + AB3)

    return X, t

#============================================================================
# To visualize the function to find the zeros
#============================================================================

def F(N, x0, start, xi, xf, step, x1, n, f):
    '''
    Compute the function to find the zeros to have only an idea of ​​where to look
    Parameters
    ----------
    N : Integer
        number of integration steps.
    x0 : float
        initial condition on position.
    start : float
        initial condition on speed.
    xi : float
        initial time of integration.
    xf : float
        final time of integration.
    step : float
        start increment
    x1 : float
        boundary condition of solution
    n : int
        number of function values ​​to calculate
    f : callable
        function to integrate, must accept vectorial input
    Returns
    -------
    xs : one dimensional array
        solution of the equation
    '''
    P = np.zeros(n)
    S = np.linspace(start, start+n*step, n)
    for j, s in enumerate(S):
        P[j] = AMB4(N, xi, xf, f, init=(x0, s))[0][-1, 0]
    return S, P-x1

#============================================================================
# Binary research to find the right solution with shooting method
#============================================================================

def SH(N, x0, start, xi, xf, step, x1, tau, f):
    '''
    Function that calculates zeros with the bisection method
    Parameters
    ----------
    N : Integer
        number of integration steps.
    x0 : float
        initial condition on position.
    start : float
        initial condition on speed.
    xi : float
        initial time of integration.
    xf : float
        final time of integration.
    step : float
        start increment
    x1 : float
        boundary condition of solution
    tau : float
        tollerance on find value
    f : callable
        function to integrate, must accept vectorial input
    Returns
    -------
    m : float
        ideal intial condition for speed
    sol : one dimensional array
        solution of the equation
    '''
    a = start
    sol = AMB4(N, xi, xf, f, init=(x0, a))
    k = sol[0][-1, 0] - x1
    while True:
        b = a + step
        sol = AMB4(N, xi, xf, f, init=(x0, b))
        D = sol[0][-1, 0] - x1
        if (k*D)<0.0:
            break
        k = D
        a = b
    while abs(a - b)>tau:
        m = (a + b)/2.0
        sol = AMB4(N, xi, xf, f, init=(x0, m))
        M = sol[0][-1, 0] - x1
        if (M*k)>0 :
            k = M
            a = m
        else :
            D = M
            b = m
    return m, sol

#============================================================================
# Main code
#============================================================================

def f(t, Y):
    x, v = Y
    x_dot = v
    v_dot = (3/2)*x**2
    return np.array([x_dot, v_dot])

xi  = 0        # left end of the interval
xf  = 1        # right end of the range
N   = 1000     # number of points
x0  = 4        # initial value at the left end
x1  = 1        # value we want assume the solution in the right extreme
tau = 1e-10    # tollerance


t, y = F(N, x0, -40, xi, xf, 0.0833, x1, 408, f)
# to visualize the zeros
plt.figure(1)
plt.title('Function to find the zeros of', fontsize=20)
plt.ylabel('x(1;s)-x(1)', fontsize=15)
plt.xlabel('s', fontsize=15)
plt.grid()
plt.plot(t, 0*t, color='red', linestyle='--')
plt.plot(t, y)

#===============================================================

v0, sol1 = SH(N, x0, -40, xi, xf, 0.1, x1, tau, f)
v1, sol2 = SH(N, x0, -10, xi, xf, 0.1, x1, tau, f)

sol, t1 = sol1
y1 , _  = sol.T 
sol, t2 = sol2
y2 , _  = sol.T

plt.figure(2)
plt.title('Two different initial conditions \n give the same boundary condition', fontsize=15)
           
print('x_1(t_f) - x1 = %e \nx_2(t_f) - x1 = %e ' %(y1[-1]-x1, y2[-1]-x1))

plt.ylabel('x(t)', fontsize=15)
plt.xlabel('t', fontsize=15)
plt.grid()
plt.plot(t1, y1, 'k', label='$x_1(t)$, $\dot{x}(t=0)$'+f'={v0:.3f}')
plt.plot(t2, y2, 'b', label='$x_2(t)$, $\dot{x}(t=0)$'+f'={v1:.3f}')
plt.legend(loc='best')
plt.show()

