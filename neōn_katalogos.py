"""
As the name might suggest this code is a list containing
several examples of implicit, explicit and symplectic
ode-solving algorithms, all of them, for simplicity are
applied to the harmonic oscillator
"""

import numpy as np
import scipy.integrate
import scipy.optimize as so
import  matplotlib.pyplot  as  plt


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


## Per alcuni integratori quanto sopra va scritto così


def F(Y, o0):
    """
    Accelerazione del sistema

    Parameters
    ----------
    Y : 1darray
        array of variables
    o0 : float
        model's parameters

    Return
    ------
    Y_ddot : 1darray
        array of force, or acceleration
    """
    x, = Y
    x_ddot = -o0*x
    Y_ddot = np.array([x_ddot])

    return Y_ddot


## Per altri ancora così



def sist(V, dt, x0, v0, o0):
    """
    Funzione per il metodo del punto medio implicito
    Ad ogni passo di integrazione va risolto il sistema

    Parameters
    ----------
    V : 1darray
        array of variables
    dt : float
        integration step
    x0, v0: float
        solution at time t
    o0 : float
        model's parameters

    Return
    ------
    list
    equations whose solution is
    the solution at time t + dt
    """
    x1, v1 = V

    R1 = v1 - v0 - dt*(-o0*(x1 + x0)/2)
    R2 = x1 - x0 - dt*(v1 + v0)/2

    return [R1, R2]


## Metodo di eluero

def eulero(num_steps, tf, f, init, args=()):
    """
    Integrator with eulwer method

    Parameters
    ----------
    num_steps : int
        number of point of solution
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


    for i in range(num_steps):
        df = f(t[i], X[i, :], *args)
        X[i + 1, :] = X[i, :] + dt*df
        t[i + 1] = t[i] + dt

    return X, t


## Metodo di eluero semi implicito (integratore simplettico)


def eulero_semi_impl(num_steps, tf, f, init, args=()):
    """
    Integrator with semi implicit eluer

    Parameters
    ----------
    num_steps : int
        number of point of solution
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


    for i in range(num_steps):
        df = f(t[i], X[i, :], *args)
        X[i + 1, 1::2] = X[i, 1::2] + dt*df[1::2]

        df = f(t[i], X[i+1, :], *args)
        X[i + 1, ::2] = X[i, ::2] + dt*df[::2]

        t[i + 1] = t[i] + dt

    return X, t


## Metodo velocity verlet (integratore simplettico)


def vel_ver(num_steps, tf, f, init, args=()):
    """
    Integrator with velocity-verlet order method

    Parameters
    ----------
    num_steps : int
        number of point of solution
    tf : float
        upper bound of integration
    f : callable
        acceleration on the system
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

    for i in range(num_steps):  #ciclo sui tempi

        acc1 = f(X[i, ::2], *args)
        X[i + 1, ::2] = X[i, ::2] + dt*X[i, 1::2] + 0.5*acc1*dt**2
        acc2 = f(X[i+1, ::2], *args)
        X[i + 1, 1::2] = X[i, 1::2] + 0.5*(acc1+acc2)*dt

        t[i + 1] = t[i] + dt

    return X, t


## Metodo runge-kutta di ordine 4


def RK4(num_steps, tf, f, init, args=()):
    """
    Integrator With Ruge-Kutta 4th order method

    Parameters
    ----------
    num_steps : int
        number of point of solution
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

    #primi passi con runge kutta
    for i in range(num_steps):
        xk1 = f(t[i], X[i, :], *args)
        xk2 = f(t[i] + dt/2, X[i, :] + xk1*dt/2, *args)
        xk3 = f(t[i] + dt/2, X[i, :] + xk2*dt/2, *args)
        xk4 = f(t[i] + dt, X[i, :] + xk3*dt, *args)
        X[i + 1, :] = X[i, :] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
        t[i + 1] = t[i] + dt

    return X, t


##Soluzione numerica con il metodo del punto medio implicito (integratore simplettico)


def implicit_mid_point(num_steps, tf, f, init, args=()):
    """
    Integrator with implicit mid point method
    We must solve a syestem of equation so do
    it with scipy.optimize

    Parameters
    ----------
    num_steps : int
        number of point of solution
    tf : float
        upper bound of integration
    f : callable
        function to integrate, or better system to solve
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

    for i in range(num_steps):
        xstart = X[i, :]
        X[i + 1, :] = so.fsolve(f, xstart, args=(dt, *X[i, :], *args))
        t[i + 1] = t[i] + dt
    return X, t


## Metodo Yoshida 4-ordine (integratore simplettico)



def Yoshida4(num_steps, tf, f, init, args=()):
    """
    Integrator with Yoshida method

    Parameters
    ----------
    num_steps : int
        number of point of solution
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
    #some funny coefficents
    l = 2**(1/3)
    w0 = -l/(2-l)
    w1 = 1/(2-l)
    #other funny coefficents
    c1 = c4 = w1/2
    c2 = c3 = (w0 + w1)/2
    d1 = d3 = w1
    d2 = w0
    #time steps
    dt = tf/num_steps
    X = np.zeros((num_steps + 1, len(init))) #matrice delle soluzioni
    t = np.zeros(num_steps + 1)              #array dei tempi

    X[0, :] = init                           #condizioni iniziali

    for i in range(num_steps):               #ciclo sui tempi
        x0 = X[i, ::2]
        v0 = X[i,1::2]

        x1 = x0 + c1*v0*dt
        v1 = v0 + d1*f(x1, *args)*dt

        x2 = x1 + c2*v1*dt
        v2 = v1 + d2*f(x2, *args)*dt

        x3 = x2 + c3*v2*dt
        v3 = v2 + d3*f(x3, *args)*dt

        X[i + 1, ::2] = x3 + c4*v3*dt
        X[i + 1,1::2] = v3
        t[i + 1] = t[i] + dt

    return X, t


## Predizione- correzzione usando il metodo di eulero e quello dei trapezi


def PC(num_steps, tf, f, init, args=()):
    """
    Integrator with predictor-corrector method
    euler and trapezoidal rule

    Parameters
    ----------
    num_steps : int
        number of point of solution
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

    for i in range(num_steps):
        #predico
        df1 = f(t[i], X[i, :], *args)
        X[i + 1, :] = X[i, :] + dt*df1
        t[i + 1] = t[i] + dt
        #corrggo
        df2 = f(t[i+1], X[i+1, :], *args)
        X[i + 1, :] = X[i, :] + 0.5*dt*(df1 + df2)

    return X, t


## Predittore correttore di ordine 4 Adamas-Bashforth-Moulton


def AMB4(num_steps, tf, f, init, args=()):
    """
    Integrator with Adams-Bashforth-Moulton
    predictor and corretor of order 4

    Parameters
    ----------
    num_steps : int
        number of point of solution
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


## Metodo del punto medio esplicito (integratore simplettico)


def mid_point(num_steps, tf, f, init, args=()):
    """
    Integrator with mid_point rule

    Parameters
    ----------
    num_steps : int
        number of point of solution
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

    for i in range(num_steps):
        xk1 = f(t[i], X[i, :], *args)
        xk2 = f(t[i], X[i, :] + xk1*dt/2, *args)
        X[i + 1, :] = X[i, :] + dt*xk2
        t[i + 1] = t[i] + dt

    return X, t


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
    num_steps =  10000

    #odeint
    ts0 = np.linspace(0, tf, num_steps + 1)
    sol = scipy.integrate.odeint(osc, init, ts0, args=(o0,), tfirst=True)
    xs0, vs0 = sol.T
    #eulero
    sol, ts1 = eulero(num_steps, tf, osc, init, args=(o0,))
    xs1, vs1 = sol.T
    #eulero semi implicito
    sol, ts2 = eulero_semi_impl(num_steps, tf, osc, init, args=(o0,))
    xs2, vs2 = sol.T
    #velocity verlet
    sol, ts3 = vel_ver(num_steps, tf, F, init, args=(o0,))
    xs3, vs3 = sol.T
    #runge kutta 4
    sol, ts4 = RK4(num_steps, tf, osc, init, args=(o0,))
    xs4, vs4 = sol.T
    #punto medio implicito
    sol, ts5 = implicit_mid_point(num_steps, tf, sist, init, args=(o0,))
    xs5, vs5 = sol.T
    #Yoshida
    sol, ts6 = Yoshida4(num_steps, tf, F, init, args=(o0,))
    xs6, vs6 = sol.T
    #Pedittore, correttore
    sol, ts7 = PC(num_steps, tf, osc, init, args=(o0,))
    xs7, vs7 = sol.T
    #Pedittore, correttore quarto ordine Adamas-Bashforth-Moulton
    sol, ts8 = AMB4(num_steps, tf, osc, init, args=(o0,))
    xs8, vs8 = sol.T
    #punto medio esplicito
    #Pedittore, correttore quarto ordine Adamas-Bashforth-Moulton
    sol, ts9 = mid_point(num_steps, tf, osc, init, args=(o0,))
    xs9, vs9 = sol.T


    ##Grafico soluzioni


    plt.figure(1)
    plt.title('Confronto soluzioni', fontsize=20)
    plt.xlabel('t', fontsize=15)
    plt.ylabel(r'$\vartheta(t)$', fontsize=15)
    plt.plot(ts0, Sol(ts0, o0, *init), 'black', label='sol analitica')
    plt.plot(ts0, xs0, 'blue', label='odeint')
    plt.plot(ts1, xs1, 'red', label='Eulero')
    plt.plot(ts2, xs2, 'green', label='Eulero semi implicito (integratore simplettico)')
    plt.plot(ts3, xs3, 'yellow', label='velocity verlet (integratore simplettico)')
    plt.plot(ts4, xs4, 'pink', label='Runge Kutta 4')
    plt.plot(ts5, xs5, 'orange', label='punto medio implicito (integratore simplettico)')
    plt.plot(ts6, xs6, 'fuchsia', label='Yoshida 4')
    plt.plot(ts7, xs7, 'violet', label='pred-corr')
    plt.plot(ts8, xs8, 'khaki', label='PCAMB4')
    plt.plot(ts9, xs9, 'navy', label='punto medio esplicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()


    ##Grafico differenze


    plt.figure(3)
    plt.suptitle('Differenza tra soluzione esatta e numerica', fontsize=20)

    plt.subplot(511)
    plt.plot(ts0, Sol(ts0, o0, *init)-xs0, 'k', label='odeint')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(512)
    plt.plot(ts1, Sol(ts1, o0, *init)-xs1, 'k', label='Eulero')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(513)
    plt.plot(ts2, Sol(ts2, o0, *init)-xs2, 'k', label='Eulero semi implicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(514)
    plt.plot(ts3, Sol(ts3, o0, *init)-xs3, 'k', label='velocity verlet (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(515)
    plt.plot(ts4, Sol(ts4, o0, *init)-xs4, 'k', label='Runge Kutta 4')
    plt.legend(loc='best')
    plt.grid()


    plt.figure(4)
    plt.suptitle('Differenza tra soluzione esatta e numerica', fontsize=20)

    plt.subplot(511)
    plt.plot(ts5, Sol(ts5, o0, *init)-xs5, 'k', label='punto medio implicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(512)
    plt.plot(ts6, Sol(ts6, o0, *init)-xs6, 'k', label='Yoshida 4(integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(513)
    plt.plot(ts7, Sol(ts7, o0, *init)-xs7, 'k', label='pred-corr')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(514)
    plt.plot(ts8, Sol(ts8, o0, *init)-xs8, 'k', label='PCAMB4')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(515)
    plt.plot(ts9, Sol(ts9, o0, *init)-xs9, 'k', label='punto medio esplicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()


    ##Grafico dell'energia


    def U(v, x):
        return (v**2 + o0*x**2)-(v[0]**2 + o0*x[0]**2)


    plt.figure(5)
    plt.suptitle('Differenza fra enerigia iniziale  ed energia al tempo t del sistema', fontsize=20)

    plt.subplot(511)
    plt.plot(ts0, U(vs0, xs0) , 'k', label='odeint')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(512)
    plt.plot(ts1, U(vs1, xs1), 'k', label='Eulero')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(513)
    plt.plot(ts2, U(vs2, xs2), 'k', label='Eulero semi implicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(514)
    plt.plot(ts3, U(vs3, xs3), 'k', label='velocity verlet (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(515)
    plt.plot(ts4, U(vs4, xs4), 'k', label='Runge Kutta 4')
    plt.legend(loc='best')
    plt.grid()


    plt.figure(6)
    plt.suptitle('Differenza fra enerigia iniziale  ed energia al tempo t del sistema', fontsize=20)

    plt.subplot(511)
    plt.plot(ts5, U(vs5, xs5), 'k', label='punto medio implicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(512)
    plt.plot(ts6, U(vs6, xs6), 'k', label='Yoshida 4 (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(513)
    plt.plot(ts7, U(vs7, xs7), 'k', label='pred-corr')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(514)
    plt.plot(ts8,  U(vs8, xs8), 'k', label='PCAMB4')
    plt.legend(loc='best')
    plt.grid()

    plt.subplot(515)
    plt.plot(ts9,  U(vs9, xs9), 'k', label='punto medio esplicito (integratore simplettico)')
    plt.legend(loc='best')
    plt.grid()


    plt.show()