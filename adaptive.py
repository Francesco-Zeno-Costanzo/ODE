import numpy as np
import matplotlib.pyplot as plt

## Runge–Kutta–Fehlberg method

def RKF45(ti, tf, f, init, tol, args=()):
    """
    Integrator with Runge-Kutta-Fehlberg method

    Parameters
    ----------
    ti : float
        lower bound of integration
    tf : float
        upper bound of integration
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
    Y : array, shape (i, len(init))
        solution of equation
    t : 1darray
        time
    H : 1darray
        array of steps
    """

    #initial time steps
    dt = abs(tf - ti)/num_steps

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
    B1 = np.array([0, 2/9, 1/12, 69/128, -17/12, 65/432]) # For y part of f
    B2 = np.array([0, 0,   1/4, -243/128, 27/4,  -5/16])
    B3 = np.array([0, 0,   0,    135/64, -27/5,  13/16])
    B4 = np.array([0, 0,   0,    0,      16/15,  4/27])
    B5 = np.array([0, 0,   0,    0,      0,      5/144])

    CH = np.array([47/450, 0, 12/25, 32/225, 1/30,  6/25]) # For new solution
    CT = np.array([-1/150, 0, 3/100, -16/75, -1/20, 6/25]) # For error computation

    iter = 0
    while time < tf:

        if time + h > tf:
            h = tf - time

        k1 = h*f(time+A[0]*h, y, *args)
        k2 = h*f(time+A[1]*h, y + B1[1]*k1, *args)
        k3 = h*f(time+A[2]*h, y + B1[2]*k1 + B2[2]*k2, *args)
        k4 = h*f(time+A[3]*h, y + B1[3]*k1 + B2[3]*k2 + B3[3]*k3, *args)
        k5 = h*f(time+A[4]*h, y + B1[4]*k1 + B2[4]*k2 + B3[4]*k3 + B4[4]*k4, *args)
        k6 = h*f(time+A[5]*h, y + B1[5]*k1 + B2[5]*k2 + B3[5]*k3 + B4[5]*k4 + B5[5]*k5, *args)

        # Truncation error, we consider the mean for all equations
        TE = abs(np.mean(CT[0]*k1 + CT[1]*k2+ CT[2]*k3+ CT[3]*k4+ CT[4]*k5+ CT[5]*k6))

        if TE < tol:
            y     = y + CH[0]*k1 + CH[1]*k2+ CH[2]*k3+ CH[3]*k4+ CH[4]*k5+ CH[5]*k6
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


def solve(t0, tf, h_t, init, tau, derivs, args=()):
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
            errarr = [yerr[i]/yscal[i] for i in range(len(yerr))]
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
        
        # Seems problematic in some cases
        if errmax > ERRCON:
            hnext = SAFETY * h * errmax**PGROW
        else:
            hnext = h*5

        hdid = h
        x += hdid
        y = ytemp

        return x, y, hdid, hnext


    #Cash-Karp Runge-Kutta step
    def rkck(y, dydx, x, h, derivs, args=()):
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
        derivs : callable
            function to integrate, must accept vectorial input
        args : tuple, optional
            extra arguments to pass to derivs

        Return
        ------
        yout : 1darray
            solution
        yerr : 1darray
            error
        """
        # For temporal part of f
        a2 = 0.2 ; a3 = 0.3 ; a4 = 0.6 ; a5 = 1.0 ; a6 = 0.875

        # For solution part of f
        b21 = 0.2           ; b31 = 3.0/40.0         ; b32 = 9.0/40.0
        b41 = 0.3           ; b42 = -0.9             ; b43 = 1.2
        b51 = -11.0/54.0    ; b52 = 2.5              ; b53 = -70.0/27.0
        b54 = 35.0/27.0     ; b61 = 1631.0/55296.0   ; b62 = 175.0/512.0
        b63 = 575.0/13824.0 ; b64 = 44275.0/110592.0 ; b65 = 253.0/4096.0

        # For update the solutions
        c1 = 37.0/378.0  ; c3 = 250.0/621.0
        c4 = 125.0/594.0 ; c6 = 512.0/1771.0

        # For the computattion of truncation error
        dc5 = -277.00/14336.0
        dc1 = c1-2825.0/27648.0  ; dc3 = c3-18575.0/48384.0
        dc4 = c4-13525.0/55296.0 ; dc6 = c6-0.25

        # step of algorithm
        y_temp = y + b21*h*dydx               # First step
        k2 = derivs(x + a2*h, y_temp, *args)  # Second step

        y_temp = y + h*(b31*dydx + b32*k2)
        k3 = derivs(x + a3*h, y_temp, *args)  # Third step

        y_temp = y + h*(b41*dydx + b42*k2 + b43*k3)
        k4 = derivs(x + a4*h, y_temp, *args)  # Fourth step

        y_temp = y + h*(b51*dydx + b52*k2 + b53*k3 + b54*k4)
        k5 = derivs(x + a5*h, y_temp, *args)  # Fifth step.

        y_temp = y + h*(b61*dydx + b62*k2 + b63*k3 + b64*k4 + b65*k5)
        k6 = derivs(x + a6*h, y_temp, *args)  # Sixth step

        # Ffinal solution, accumulate increments with proper weights
        yout = y + h*(c1*dydx + c3*k3 + c4*k4 + c6*k6)
        
        # Estimate error as difference between fourth and fifth order methods
        yerr = h*(dc1*dydx + dc3*k3 + dc4*k4 + dc5*k5 + dc6*k6)

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
        #return -2/(t**2 -2)
        #return (1-t)*np.exp(t)
        return p[2]/np.sqrt(p[0])*np.sin(np.sqrt(p[0])*t) + p[1]*np.cos(np.sqrt(p[0])*t)

    def eq(t, Y, o0):
        """ Equation to solve
        """
        x, v = Y
        #-----------------
        #x_dot = v
        #v_dot = 2*v - x
        #-----------------
        x_dot = v
        v_dot = - o0 * x
        #-----------------
        Y_dot = np.array([x_dot, v_dot])

        #x, = Y
        #x_dot = x**2 * t
        #Y_dot = np.array([x_dot])        

        return Y_dot

    # Parameter of equations
    o0 = 9
    # Initial conditions
    v0 = 0
    x0 = 1
    init = np.array([x0 , v0]) #x(0), x(0)'
    # bound of interval
    ti = 0
    tf = 10

    sol, ts0, hs0 = RKF45(ti, tf, eq, init, 1e-12, args=(o0,))
    xs0, _ = sol.T
    sol, ts1, hs1 = solve(ti, tf, 0.05, init, 1e-12, eq, args=(o0, ))
    xs1, _ = sol.T

    # Solutions
    ts3 = np.linspace(ti, tf, 1e4)
    par = np.array([o0, *init])
    plt.figure(1)
    plt.title("Soluzione con algoritmo adattivo", fontsize=15)
    plt.grid()
    plt.plot(ts0, xs0, label='RKF45')
    plt.plot(ts1, xs1, label='Cash-Karp')
    plt.plot(ts3, Sol(ts3, par), 'k', label='exact')
    plt.legend(loc='best')

    # Global Error
    plt.figure(2)
    plt.suptitle('Differenza tra soluzione esatta e numerica', fontsize=20)

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
    plt.suptitle('Paso di integrazione', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, hs0, 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    plt.plot(ts1, hs1, 'k', label='Cash-Karp')
    plt.legend(loc='best')
    plt.grid()

    plt.show()