import numpy as np
import matplotlib.pyplot as plt


## Soluzione analitica


def Sol(t, o0, x0, v0):
    """Analitic solutions
    """
    return v0/np.sqrt(o0) * np.sin(np.sqrt(o0)*t) + x0*np.cos(np.sqrt(o0)*t)


##Equazioni da risolvere


def osc(t, Y, o0):
    """
    equation to solve

    Parameters
    ----------
    t : float
        time
    Y : 1darray
        array of variables
    o0 : float
        model's parameters

    Return
    ------
    Y_dot : 1darray
        array of equations
    """
    theta, omega = Y

    theta_dot = omega
    omega_dot = -o0 * theta

    Y_dot = np.array([theta_dot, omega_dot])

    return Y_dot


## Runge–Kutta–Fehlberg method


def RKF45(num_steps, tf, f, init, tau, args=()):
    """
    Integrator with Runge–Kutta–Fehlberg method

    Parameters
    ----------
    num_steps : int
        max number of point of solution
    tf : float
        upper bound of integration
    f : callable
        function to integrate, must accept vectorial input
    init : 1darray
        array of initial condition
    tau : float
        required accuracy
    args : tuple, optional
        extra arguments to pass to f

    Return
    ------
    X : array, shape (i, len(init))
        solution of equation
    t : 1darray
        time
    h : 1darray
        array of steps
    """

    #initial time steps
    dt = tf/num_steps

    X = np.zeros((num_steps + 1, len(init))) #matrice delle soluzioni
    t = np.zeros(num_steps + 1)              #array dei tempi
    h = np.zeros(num_steps + 1)              #array dei passi

    X[0, :] = init                           #condizioni iniziali
    h[0] = dt
    i = 0

    A  = np.array([0, 2/9, 1/3,  3/4,    1,      5/6])
    B1 = np.array([0, 2/9, 1/12, 69/128, -17/12, 65/432])
    B2 = np.array([0, 0,   1/4, -243/128, 27/4,  -5/16])
    B3 = np.array([0, 0,   0,    135/64, -27/5,  13/16])
    B4 = np.array([0, 0,   0,    0,      16/15,  4/27])
    B5 = np.array([0, 0,   0,    0,      0,      5/144])

    CH = np.array([47/450, 0, 12/25, 32/225, 1/30,  6/25])
    CT = np.array([-1/150, 0, 3/100, -16/75, -1/20, 6/25])

    iter = 0
    while t[i] < tf:

        t[i+1] = t[i] + h[i]

        k1 = h[i]*f(t[i]+A[0]*h[i], X[i, :], *args)

        k2 = h[i]*f(t[i]+A[1]*h[i], X[i, :] + B1[1]*k1, *args)

        k3 = h[i]*f(t[i]+A[2]*h[i], X[i, :] + B1[2]*k1 + B2[2]*k2, *args)

        k4 = h[i]*f(t[i]+A[3]*h[i], X[i, :] + B1[3]*k1 + B2[3]*k2 + B3[3]*k3, *args)

        k5 = h[i]*f(t[i]+A[4]*h[i], X[i, :] + B1[4]*k1 + B2[4]*k2 + B3[4]*k3 + B4[4]*k4, *args)

        k6 = h[i]*f(t[i]+A[5]*h[i], X[i, :] + B1[5]*k1 + B2[5]*k2 + B3[5]*k3 + B4[5]*k4 + B5[5]*k5, *args)


        X[i+1, :] = X[i, :] + CH[0]*k1 + CH[1]*k2+ CH[2]*k3+ CH[3]*k4+ CH[4]*k5+ CH[5]*k6

        TE = abs(np.mean(CT[0]*k1 + CT[1]*k2+ CT[2]*k3+ CT[3]*k4+ CT[4]*k5+ CT[5]*k6))

        hn = 0.9*h[i]*(tau/TE)**(1/5)
        if TE < tau:
            i += 1
        h[i] = hn
        iter += 1


    X = X[:i, :]
    t = t[:i]
    h = h[:i]
    print(f'numero di terazioni = {iter}')

    return X, t, h

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

            htemp = SAFETY*h*pow(errmax, PSHRINK)#reducing step size, try again
            if h >= 0.0:
                h = max(htemp, 0.1*h)
            else:
                h = min(htemp, 0.1*h)

            xnew = x+h
            if xnew == x:
                print("Stepsize underflow in rkqs")

        if errmax > ERRCON:
            hnext = SAFETY*h*pow(errmax, PGROW)
        else:
            hnext = 5.0*h

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

        #tantissimi numeri
        #per il passo temporale
        a2 = 0.2 ; a3 = 0.3 ; a4 = 0.6 ; a5 = 1.0 ; a6 = 0.875

        #per le chiamate successive
        b21 = 0.2           ; b31 = 3.0/40.0         ; b32 = 9.0/40.0
        b41 = 0.3           ; b42 = -0.9             ; b43 = 1.2
        b51 = -11.0/54.0    ; b52 = 2.5              ; b53 = -70.0/27.0
        b54 = 35.0/27.0     ; b61 = 1631.0/55296.0   ; b62 = 175.0/512.0
        b63 = 575.0/13824.0 ; b64 = 44275.0/110592.0 ; b65 = 253.0/4096.0

        #per aggiornamento soluzione
        c1 = 37.0/378.0  ; c3 = 250.0/621.0
        c4 = 125.0/594.0 ; c6 = 512.0/1771.0

        #per calcolo dell'errore
        dc5 = -277.00/14336.0
        dc1 = c1-2825.0/27648.0  ; dc3 = c3-18575.0/48384.0
        dc4 = c4-13525.0/55296.0 ; dc6 = c6-0.25

        #passi dell'algoritmo
        y_temp = y + b21*h*dydx                #first step
        ak2 = derivs(x + a2*h, y_temp, *args)  #second step

        y_temp = y + h*(b31*dydx + b32*ak2)
        ak3 = derivs(x + a3*h, y_temp, *args)  #third step

        y_temp = y + h*(b41*dydx + b42*ak2 + b43*ak3)
        ak4 = derivs(x + a4*h, y_temp, *args)  #Fourth step

        y_temp = y + h*(b51*dydx + b52*ak2 + b53*ak3 + b54*ak4)
        ak5 = derivs(x + a5*h, y_temp, *args)  #Fifth step.

        y_temp = y + h*(b61*dydx + b62*ak2 + b63*ak3 + b64*ak4 + b65*ak5)
        ak6 = derivs(x + a6*h, y_temp, *args)  #Sixth step

        #final solution, accumulate increments with proper weights
        yout = y + h*(c1*dydx + c3*ak3 + c4*ak4 + c6*ak6)
        #errore
        #Estimate error as difference between fourth and fifth order methods
        yerr = h*(dc1*dydx + dc3*ak3 + dc4*ak4 + dc5*ak5 + dc6*ak6)

        return yout, yerr

    # Integrazione

    TINY = 1e-30
    Y = []                              #to store solution
    t = []                              #to store time
    H = []                              #to store integration steps
    x = t0                              #initial point of integration
    h = h_t*(tf - t0)/abs(tf - t0)      #initial step (tf - t0)/abs(tf - t0) cab be +-1
    y = init                            #initial condition

    Y.append(y)                         #store the first point
    t.append(x)                         #store the first time
    H.append(h)                         #store the first step
    iter = 0                            #to count iterations

    while x < tf :

        dydx = derivs(x, y, *args)

        yscal = abs(y) + abs(h*dydx) + TINY
        if (x+h-tf)*(x+h-t0) > 0 : h = tf-x #If stepsize can overshoot, decrease

        x, ytemp, hdid, hnext = rkqs(y, dydx, x, h, tau, yscal, derivs, args)
        H.append(hdid)     #store the steps used
        y = ytemp          #update the solution
        h = hdid#next      #update step

        Y.append(y)        #store solution
        t.append(x)        #store time
        iter += 1          #update iteration

    print(f'numero di terazioni = {iter}')

    return np.array(Y), np.array(t), np.array(H)

## Risoluzione

if __name__ == '__main__':

    #parametri simulazione
    o0 = 9
    #condizione iniziali
    v0 = 0
    x0 = 1
    init = np.array([x0 , v0]) #x(0), x(0)'
    #estremi di integrazione
    ti = 0
    tf = 10
    #numero di punti
    num_steps = int(1e5)

    sol, ts0, hs0 = RKF45(num_steps, tf, osc, init, 1e-12, args=(o0,))
    xs0, vs0 = sol.T
    sol, ts1, hs1 = solve(ti, tf, 0.05, init, 1e-12, osc, args=(o0, ))
    xs1, vs1 = sol.T

## Grafico soluzioni

    plt.figure(1)
    plt.title("Soluzione con algoritmo adattivo", fontsize=15)
    plt.grid()
    plt.plot(ts0, xs0, label='RKF45')
    plt.plot(ts1, xs1, label='Cash-Karp')
    plt.legend(loc='best')

## Grafico differenze

    plt.figure(2)
    plt.suptitle('Differenza tra soluzione esatta e numerica', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, Sol(ts0, o0, *init)-xs0, 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    plt.plot(ts1, Sol(ts1, o0, *init)-xs1, 'k', label='Cash-Karp')
    plt.legend(loc='best')
    plt.grid()

## Grafico dell'energia


    def U(v, x):
        """
        Energy of system

        Parameters
        ----------
        v : 1darray
            velocity
        x : 1darray
            position

        Return
        ------
        E(t) - E(t=0)
        """
        return (v**2 + o0*x**2)-(v[0]**2 + o0*x[0]**2)

    plt.figure(3)
    plt.suptitle('Differenza fra enerigia iniziale  ed energia al tempo t del sistema', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, U(vs0, xs0), 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    plt.plot(ts1, U(vs1, xs1), 'k', label='Cash-Karp')
    plt.legend(loc='best')
    plt.grid()

## Grafico andamento passo

    plt.figure(4)
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