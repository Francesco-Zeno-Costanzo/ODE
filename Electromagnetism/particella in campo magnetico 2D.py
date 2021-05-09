import numpy as np
import scipy.integrate
import  matplotlib.pyplot  as  plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation


def f(y, t, b, c):
    x, u, z, v = y        #vettore con le incognite
    dydt = [u, b*t*v, v, -b*t*u] #vettore parte di destra
    return dydt


b=1 #vaore dei coefficienti
c=1
x0 = [0, 10, 0, 10]     #condizioni iniziali;in ordine: x(0), x(0)', y(0), y(0)',
t = np.linspace(-10, 10, 10000)    #tempo su cui vuoi integrare (in sintesi i passetti)
sol = scipy.integrate.odeint(f, x0, t, args=(b, c))

x = sol[:, 0]
y = sol[:, 2]

#grafico 3D
fig = plt.figure(1)
ax = fig.gca(projection='3d')
ax.plot(sol[:, 0], sol[:, 2], t) # plot in ordine: x, y, z.
ax.set_title('Carica in campo magnetico lineare nel tempo')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('Tempo')


fig = plt.figure(2)

plt.xlim(np.min(x),np.max(x))
plt.ylim(np.min(y),np.max(y))

line, = plt.plot([],[],'b-')
dot, = ax.plot([], [], 'ro')

def animate(i):
    dot.set_data(x[i],y[i])
    line.set_data(x[:i],y[:i])
    return  dot, line


myAnimation = animation.FuncAnimation(fig, animate, frames=range(0, 10000, 15),interval=10, blit=True, repeat=True)

plt.grid()
plt.title('Carica in campo magnetico lineare nel tempo')
plt.xlabel('X(t)')
plt.ylabel('Y(t)')

#anim.save('pinem.mp4', fps=30, extra_args=['-vcodec', 'libx264']) #salva sulla cartella corrente l'animazione
plt.show()