import numpy as np
import matplotlib.pyplot as plt

## Runge–Kutta–Fehlberg method

def RKF45(ti, tf, h_t, f, init, tol, args=()):
    """
    Integrator with Runge-Kutta-Fehlberg method

    Parameters
    ----------
    ti : float
        lower bound of integration
    tf : float
        upper bound of integration
    h_t : float
        guess for integration steps
    f : callable
        function to integrate, must accept vectorial input
    init : 1darray
        array of initial condition
    tol: float
        required accuracy
    args : tuple, optional
        extra arguments to pass to f

    Return
    ------
    Y : array, shape (num of steps done, len(init))
        solution of equation
    t : 1darray
        time
    H : 1darray
        array of steps
    """

    #initial time steps
    dt = h_t

    Y = []         # To store solutions
    t = []         # Array of time
    H = []         # Array to store size step during integration

    Y.append(init) # Set initial condition
    t.append(ti)
    H.append(dt)
    y    = init
    time = ti
    i    = 0
    h    = dt

    # Coefficients for evolution
    A  = np.array([0, 2/9, 1/3,  3/4,    1,      5/6])    # For t part of f
    # For y part of f
    B = np.array([
        [0,       0,       0,      0,     0,     0],
        [2/9,     0,       0,      0,     0,     0],
        [1/12,    1/4,     0,      0,     0,     0],
        [69/128, -243/128, 135/64, 0,     0,     0],
        [-17/12,  27/4,   -27/5,   16/15, 0,     0],
        [65/432, -5/16,    13/16,  4/27,  5/144, 0]
    ])

    CH = np.array([47/450, 0, 12/25, 32/225, 1/30,  6/25]) # For new solution
    CT = np.array([-1/150, 0, 3/100, -16/75, -1/20, 6/25]) # For error computation

    iter = 0
    while time < tf:

        if time + h > tf:
            h = tf - time

        k = np.zeros((6, len(y)))
        for j in range(6):
            dy = sum(B[j, l] * k[l] for l in range(j))
            k[j] = h * f(time + A[j] * h, y + dy, *args)

        # Truncation error, we consider the mean for all equations
        TE = abs(np.mean(np.dot(CT, k)))

        if TE < tol:
            y     = y + np.dot(CH, k)
            time += h
            Y.append(y)
            t.append(time)
            H.append(h)
            iter += 1

        # New step size
        h = 0.9*h*(tol/TE)**(1/5)

    Y = np.array(Y)
    t = np.array(t)
    H = np.array(H)
    print(f'Number of iterartions = {iter}')

    return Y, t, H

## Cash-Karp-Runga-Kutta Method from Numerical recipes


def solve(t0, tf, h_t, derivs, init, tau, args=()):
    """
    Integrator Cash-Karp-Runga-Kutta Method

    Parameters
    ----------
    t0 : float
        lower bound of integration
    tf : float
        upper bound of integration
    h_t : float
        guess for integration steps
    init : 1darray
        array of initial condition
    tau : float
        required accuracy
    derivs : callable
        function to integrate, must accept vectorial input
    args : tuple, optional
        extra arguments to pass to derivs

    Return
    ------
    Y : array, shape (iter, len(init))
        solution of equation
    t : 1darray
        time
    h : 1darray
        array of steps
    """

    #Useful function

    #Adaptive stepsize
    def rkqs(y, dydx, x, htry, eps, yscal, derivs, args=()):
        """
        function to controll Adaptive stepsize

        Parameters
        ----------
        y : 1darray
            solution at point x
        dydx : 1darray
            array for RHS of equation
        x : float
            actual point
        htry : float
            guess for integration steps
        eps : float
            required accuracy
        yscal : 1darray
            array to controll the accuracy
        derivs : callable
            function to integrate, must accept vectorial input
        args : tuple, optional
            extra arguments to pass to derivs

        Return
        ------
        x : float
            new point fo solution
        y : 1darray
            solution at point x
        hdid : float
            value of step used
        hnext : float
            prediction of next stepsize
        """

        # numeri importanti
        SAFETY  = 0.9
        PGROW   = -0.2
        PSHRINK = -0.25
        ERRCON  = 1.89e-4 #(5/SAFETY)**(1/PGROW)

        h = htry

        while True:
            ytemp, yerr = rkck(y, dydx, x, h, derivs, args)
            errarr = [abs(yerr[i]/yscal[i]) for i in range(len(yerr))]
            errmax = np.max(errarr)
            errmax /= eps
            if errmax <= 1.0: #Step was succesful, time for next step!
                break

            htemp = SAFETY * h * errmax**PSHRINK #reducing step size, try again
            if h >= 0.0:
                h = max(htemp, 0.1*h)
            else:
                h = min(htemp, 0.1*h)

            xnew = x+h
            if xnew == x:
                print("Stepsize underflow in rkqs")
                continue
        
        if errmax > ERRCON:
            hnext = SAFETY * h * errmax**PGROW
        else:
            hnext = h*5

        hdid = h
        x += hdid
        y = ytemp

        return x, y, hdid, hnext


    #Cash-Karp Runge-Kutta step
    def rkck(y, dydx, x, h, f, args=()):
        """
        function for integration of quation

        Parameters
        ----------
        y : 1darray
            solution at point x
        dydx : 1darray
            array for RHS of equation
        x : float
            actual point
        h : float
            integration step
        f : callable
            function to integrate, must accept vectorial input
        args : tuple, optional
            extra arguments to pass to f

        Return
        ------
        yout : 1darray
            solution
        yerr : 1darray
            error
        """
        # For temporal part of f
        A = np.array([0, 0.2, 0.3, 0.6, 1.0, 0.875])

        # For solution part of f
        B = np.array([
            [0,          0,       0,         0,            0,        0],
            [0.2,        0,       0,         0,            0,        0],
            [3/40,       9/40,    0,         0,            0,        0],
            [0.3,       -0.9,     1.2,       0,            0,        0],
            [-11/54,     2.5,    -70/27,     35/27,        0,        0],
            [1631/55296, 175/512, 575/13824, 44275/110592, 253/4096, 0]
        ])

        # For update the solutions ad truncation error
        C = np.array([37/378,     0, 250/621,     125/594,     0,         512/1771])
        D = np.array([2825/27648, 0, 18575/48384, 13525/55296, 277/14336, 0.25])
        DC = C - D

        # Step of algorithm
        k    = np.zeros((6, len(y)))
        k[0] = dydx

        for j in range(1, 6):
            y_temp = y + h * sum(B[j, l] * k[l] for l in range(j))
            k[j] = f(x + A[j] * h, y_temp, *args)
        
        # Final solution, accumulate increments with proper weights
        yout = y + h * np.dot(C, k)

        # Estimate error as difference between fourth and fifth order methods
        yerr = h * np.dot(DC, k)

        return yout, yerr

    # Integration

    TINY = 1e-30                      # To avoid division by zero
    Y = []                            # To store solution
    t = []                            # To store time
    H = []                            # To store integration steps
    x = t0                            # Initial point of integration
    h = h_t*(tf - t0)/abs(tf - t0)    # Initial step (tf - t0)/abs(tf - t0) cab be +-1
    y = init                          # Initial condition

    Y.append(y)                       # Store the first point
    t.append(x)                       # Store the first time
    H.append(h)                       # Store the first step
    iter = 0                          # To count iterations

    while x < tf :

        dydx = derivs(x, y, *args)

        yscal = abs(y) + abs(h*dydx) + TINY
        if (x+h-tf)*(x+h-t0) > 0 : h = tf-x # If stepsize can overshoot, decrease

        x, ytemp, hdid, hnext = rkqs(y, dydx, x, h, tau, yscal, derivs, args)
        H.append(hdid)     # Store the steps used
        y = ytemp          # Update the solution
        h = hnext          # Update step (sometimes hdid?? why??)

        Y.append(y)        # Store solution
        t.append(x)        # Store time
        iter += 1          # Update iteration

    print(f'Number of iterartions = {iter}')

    return np.array(Y), np.array(t), np.array(H)

## Risoluzione

if __name__ == '__main__':

    def Sol(t, p):
        """Analitic solutions
        """
        #return (1-t)*np.exp(t)
        #return p[2]/np.sqrt(p[0])*np.sin(np.sqrt(p[0])*t) + p[1]*np.cos(np.sqrt(p[0])*t)
        #return -2/(t**2 -2)
        return np.sin(2*t)*np.cos(t)/(1+(t-3)**2)


    def eq(t, Y, o0):
        """ Equation to solve
        """
        #x, v = Y
        #-----------------
        #x_dot = v
        #v_dot = 2*v - x
        #-----------------
        #x_dot = v
        #v_dot = - o0 * x
        #-----------------
        #Y_dot = np.array([x_dot, v_dot])

        x = Y
        #-----------------
        #x_dot = x**2 * t
        #-----------------
        d  = 1 + (t-3)**2
        a1 = -d*np.sin(t)*np.sin(2*t)
        a2 =  d*np.cos(t)*np.cos(2*t)*2
        a3 =   -np.cos(t)*np.sin(2*t)*2*(t-3)
        x_dot = (a1+a2+a3)/d**2
        #-----------------
        Y_dot = np.array([x_dot])        

        return Y_dot

    # Parameter of equations
    o0 = 9
    # Initial conditions
    v0 = 0
    x0 = 1
    init = np.array([0 ])#, v0]) #x(0), x(0)'
    # Bound of interval
    ti = 0
    tf = 10

    sol, ts0, hs0 = RKF45(ti, tf, 0.01, eq, init, 1e-12, args=(o0,))
    xs0, = sol.T
    sol, ts1, hs1 = solve(ti, tf, 0.01, eq, init, 1e-12, args=(o0,))
    xs1, = sol.T

    # Solutions
    ts3 = np.linspace(ti, tf, int(1e4))
    par = np.array([o0, *init])
    plt.figure(1)
    plt.title("Solutions Comparison", fontsize=15)
    plt.grid()
    plt.plot(ts0, xs0, label='RKF45')
    plt.plot(ts1, xs1, label='Cash-Karp')
    plt.plot(ts3, Sol(ts3, par), 'k', label='analytical')
    plt.legend(loc='best')

    # Global Error
    plt.figure(2)
    plt.suptitle('analytical - numerical', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, Sol(ts0, par)-xs0, 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    plt.plot(ts1, Sol(ts1, par)-xs1, 'k', label='Cash-Karp')
    plt.legend(loc='best')
    plt.grid()

    # Size step evolution
    plt.figure(3)
    plt.suptitle('Integration steps', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, hs0, 'k', label='RKF45')
    plt.legend(loc='best')
    plt.yscale('log')
    plt.grid()

    plt.subplot(122)
    plt.plot(ts1, hs1, 'k', label='Cash-Karp')
    plt.legend(loc='best')
    plt.yscale('log')
    plt.grid()

    plt.show()