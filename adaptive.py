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

## Cash-Karp Runga-Kutta Method


## Risoluzione

if __name__ == '__main__':
    #parametri simulazione
    o0 = 9
    #condizione iniziali
    v0 = 0
    x0 = 1
    init = np.array([x0 , v0]) #x(0), x(0)'
    #estremo di integrazione
    tf = 10
    #numero di punti
    num_steps = int(1e5)

    sol, ts0, hs0 = RKF45(num_steps, tf, osc, init, 1e-13, args=(o0,))
    xs0, vs0 = sol.T

## Grafico soluzioni

    plt.figure(1)
    plt.title("Soluzione con alhoritmo adattivo", fontsize=15)
    plt.grid()
    plt.plot(ts0, xs0, label='RKF45')
    plt.legend(loc='best')

## Grafico differenze

    plt.figure(2)
    plt.suptitle('Differenza tra soluzione esatta e numerica', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, Sol(ts0, o0, *init)-xs0, 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    #plt.plot(ts1, Sol(ts1, o0, *init)-xs1, 'k', label='Eulero')
    #plt.legend(loc='best')
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
    plt.plot(ts0, U(vs0, xs0) , 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    #plt.plot(ts1, U(vs1, xs1), 'k', label='Eulero')
    #plt.legend(loc='best')
    plt.grid()

## Grafico andamento passo

    plt.figure(4)
    plt.suptitle('Paso di integrazione', fontsize=20)

    plt.subplot(121)
    plt.plot(ts0, hs0 , 'k', label='RKF45')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(122)
    #plt.plot(ts1, U(vs1, xs1), 'k', label='Eulero')
    #plt.legend(loc='best')
    plt.grid()

    plt.show()