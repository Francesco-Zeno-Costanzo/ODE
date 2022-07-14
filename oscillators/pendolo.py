import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

##CAVEAT:
#Il metodo periodo della classe PendoloNumerico
#è stato scritto ipotizzando che le condizioni inizali
#siano theta(0) = theta_0 e omega(0) = 0
#quindi funziona solo in quel caso
##

l = 1
g = 9.81
o0 = g/l


class PendoloAnalitico:
    """
    Classe per il calcolo del periodo tramite integrale ellittico.
    Chiaramente l'integrale viene calcolato numericamente quindi
    c'è comunque un certo errore che si aggiunge ai problemi che
    la funzione integranda da intorno all'estremo superiore
    """

    def __init__(self, theta0):
        """condizione inziale per l'angolo
        """
        self.theta0 = theta0


    def eq_P(self, x):
        """
        Funzione integranda dell'integrale ellittico
        """
        dtdx = 1/(np.sqrt(2*(np.cos(x)-np.cos(self.theta0))))

        return dtdx


    def RK4(self, f):
        """
        Integratore con runge kutta 4 che
        fondamentalmente è cavalieri simpson
        Parametri
        ---------
        f : callable
            funzione da integrare

        Return
        ---------
        x : float
            valore dell'integrale
        """
        num_steps = 100000
        dt = self.theta0/num_steps

        x = [0.]
        t = [0.]
        #self.theta0-2*dt perchè sennò non funziona benissimo
        #dati i problemi di divergenza dell'integranda
        #quindi ci si ferma un po' prima
        while t[-1] < self.theta0-2*dt:
            xk1 = f(t[-1])
            xk2 = f(t[-1] + dt/2)
            xk3 = f(t[-1] + dt/2)
            xk4 = f(t[-1] + dt)

            x_n = x[-1] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
            t_n = t[-1] + dt

            x.append(x_n)
            t.append(t_n)

        return np.array(x[-1])


    def periodo(self):
        """
        Periodo del pendolo corretto
        """
        T0 = 2*np.pi/np.sqrt(o0) #valore per le piccole oscillazioni
        return 2*T0*self.RK4(self.eq_P)/np.pi



class PendoloNumerico:
    """
    Classe per la risoluzione numerica dell'equazione
    del moto per il pendolo semplice
    """

    def __init__(self):
        pass

    def eq_moto(self, y):
        """
        sistema da risolvere
        """
        x, v = y
        Y_dot = np.array([v, -o0*np.sin(x)])

        return Y_dot



    def RK4(self, f, y0, tf, num_steps):
        """
        Calcola la soluzione usando runge_kutta4.

        Parametri
        ---------
        f : callable
            funzione da integrare.
        y0 : 1d array
            array delle condizioni iniziali.
        tf : float
            fino a dove integrare nel tmepo.
        num_steps : float
            numero di punti
        Returns
        ---------
        X : 2d array
            matrice contenente l'angolo e la velocità angolare
        t : 1d array
            tempo
        """

        dt = tf/num_steps

        X = np.zeros((num_steps + 1, 2))
        t = np.zeros(num_steps + 1)

        X[0, :] = y0


        for i in range(num_steps):
            xk1 = f(X[i, :])
            xk2 = f(X[i, :] + xk1*dt/2)
            xk3 = f(X[i, :] + xk2*dt/2)
            xk4 = f(X[i, :] + xk3*dt)
            X[i + 1, :] = X[i, :] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
            t[i + 1] = t[i] + dt

        return X, t

    def solve(self, y0, tf=10, num_steps=10000):
        """
        Calcola la soluzione usando runge_kutta4.

        Parametri
        ---------
        f : callable
            funzione da integrare.
        y0 : 1d array
            array delle condizioni iniziali.
        tf : float, optional, default=10
            fino a dove integrare nel tmepo.
        num_steps : float, optional, default=10000
            numero di punti
        Returns
        ---------
        X : 2d array
            matrice contenente l'angolo e la velocità angolare
        t : 1d array
            tempo
        """
        Y, t = self.RK4(self.eq_moto, y0, tf=tf, num_steps=num_steps)
        return Y, t

    def periodo(self, x, t):
        """
        Data la soluzione numerica trova il periodo
        CAVEAT:
        Si è assunto, per semplicità, che la simulazione sia
        lanciata con condizione iniziale sulla velocità nulla
        quindi per calcolare il periodo si trova l'indice del
        secondo massimo di quella che è theta(t), poichè per
        costruzione il primo massimo e all'inizio dei tempi.
        Per condizioni inziali diverse questa funzione non
        restitusce il valore corretto del periodo.
        Mi perdonerte la pigrizia
        """
        t1=0
        #trovo il massimo scorrendo l'array
        for i in range(1,len(x)-1):
            if x[i]>x[i+1] and x[i]>x[i-1]:
                t1 = i
                break
        #moltiplico per il passo temporale per avere il periodo in secondi
        T = t1*t[1]
        return T

    def plot_anim(self, theta, t):
        """
        Funzione che fa il plot di theta(t)
        e l'animazione del pendolo
        Parametri
        ---------
        theta : 1d array
            angolo in funzione del tempo
        t : 1d array
            tempo
        """

        dt = t[1]
        x = l*np.sin(theta)
        y = -l*np.cos(theta)

        fig = plt.figure(1, figsize=(10, 6))
        plt.suptitle('Pendolo semplice')
        ax = fig.add_subplot(121)
        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        plt.xlim(-2, 2)
        plt.ylim(-2, 2)
        plt.gca().set_aspect('equal', adjustable='box')

        xf, yf = [0,x[0]],[0,y[0]]

        line1, = plt.plot(xf, yf, linestyle='-', marker='o',color='k')

        plt.grid()

        def animate(i):

            xf[1] = x[i]
            yf[1] = y[i]
            line1.set_data(xf, yf)
            time_text.set_text(time_template % (i*dt))

            return line1, time_text

        anim=animation.FuncAnimation(fig, animate, frames=range(0, len(t), 5), interval=1, blit=True, repeat=True)

        plt.subplot(122)
        plt.ylabel(r'$\theta$(t) [rad]')
        plt.xlabel('t [s]')
        plt.plot(t, theta)
        plt.grid()

        plt.show()


theta0 = np.pi/2
P0 = PendoloNumerico()
P1 = PendoloAnalitico(theta0)
Y, t = P0.solve([theta0, 0])
xs = Y[:,0]
p0 = P0.periodo(xs, t)
P0.plot_anim(xs, t)
p1 = P1.periodo()
print(p0, p1)
"""

T0 = []
T1 = []
theta = np.linspace(0.05, 3, 60)
a = np.ones(len(theta))
for theta0 in theta:
    P0 = PendoloNumerico()
    P1 = PendoloAnalitico(theta0)
    Y, t = P0.solve([theta0, 0])
    xs = Y[:,0]
    p0 = P0.periodo(xs, t)
    p1 = P1.periodo()
    T0.append(p0)
    T1.append(p1)
##
T0 = np.array(T0)
T1 = np.array(T1)

plt.figure(1)

plt.subplot(211)
plt.title("Periodo del pendolo in funzione dell'ampiezza", fontsize=12)
plt.plot(theta, a*2*np.pi/np.sqrt(o0), label='piccole oscilalzioni')
plt.plot(theta, T0, label = 'soluzione numerica')
plt.plot(theta, T1, label='integrale ellittico')
plt.legend(loc='best')
plt.ylabel('T [s]')
plt.grid()

plt.subplot(212)
plt.title("Differenza fra simulazione e formula analitica", fontsize=12)
plt.plot(theta, T0-T1)
plt.xlabel(r'$\theta$ [rad]')
plt.ylabel('T [s]')
plt.grid()
plt.show()
"""