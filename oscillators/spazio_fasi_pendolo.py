import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt

#numero di condizioni inziali per posizione e volocità
N = 15
M = 15

#o0 = l/g e tf = tempo di integrazione
o0 = 1
tf = 10

#range dei valori inziali
x_min = -10
x_max = -x_min
v_min = -3
v_max = -v_min

#array contenente i valori inziali
x_0 = np.linspace(x_min, x_max, N)
v_0 = np.linspace(v_min, v_max, M)

#numero di passi e grandezza del passo
num_steps = 10000
dt = tf/num_steps

def eq(x, v):
    '''sistema da risolvere
    '''

    x_dot = v
    v_dot = -o0*np.sin(x)

    return x_dot, v_dot

def F(V, x0, v0):
    '''sistema da risolvere ad ogni passo dell'integrazione
    '''

    v1, x1 = V
    dx, dv = eq((x1 + x0)/2, (v1 + v0)/2)
    R1 = v1 - v0 - dt*dv
    R2 = x1 - x0 - dt*dx

    return [R1, R2]

def sol(x0, v0):
    '''
    resistiuisce gli array soluzione
    dell'equazione del moto
    '''

    x = np.zeros(num_steps+1)
    v = np.zeros(num_steps+1)

    x[0], v[0] = x0, v0

    for i in range(num_steps):
        start = (v[i], x[i])
        v[i + 1], x[i + 1] = so.fsolve(F , start, args=(x[i], v[i]), xtol=1e-06)

    return x, v

plt.figure(1)
plt.title("Spazio delle fasi pendolo", fontsize=15)
plt.grid()
plt.xlim(x_min, x_max)

#risolvo per ogni condizione inizale per
#avere l'intero spazio delle fasi
for x0 in x_0:
    for v0 in v_0:
        theta, ptheta = sol(x0, v0)
        plt.plot(theta, ptheta, 'b')

plt.show()