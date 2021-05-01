#imoprtiamo una balla di roba
import numpy as np
import scipy.integrate
import  matplotlib.pyplot  as  plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation

'''
Bisogna portare le equazioni, se di ordine superiore, a equazioni al primo ordine affinchè siano risolte
ESEMPIO:
Sistema da dover risolvere:
x(t)" - b*y(t)' = 0
y(t)" + b*x(t)' = 0
Abbassiamo di ordine e scriviamo le in forma normale:
u(t) = x(t)'
v(t) = y(t)'
u(t)' = b*v(t)
v(t)' =-b*u(t)
che poi saranno messe in un unico vettore
'''

T=1000 #numero di passi
def f(q, t, B1, B2, B3, c, E1, E2, E3):
    x, u, y, v, z, w = q        #vettore con le incognite
    dqdt = [u, c*(B3*v-w*B2)+c*E1, v, c*(w*B1-u*B3)+c*E2, w, c*(u*B2-v*B1)+c*E3] #vettore parte di destra
    return dqdt

B1=0 #vaore dei coefficienti
B2=0
B3=1
c=1
E1=1
E2=0
E3=0
x0 = [0, 10, 0, 10, 0, 10]     #condizioni iniziali;in ordine: x(0), x(0)', y(0), y(0)', z(0), z(0)'
t = np.linspace(0, 20, T)    #tempo su cui vuoi integrare (in sintesi i passetti)
sol = scipy.integrate.odeint(f, x0, t, args=(B1, B2, B3, c,  E1, E2, E3))

x = sol[:, 0]
y = sol[:, 2]
z = sol[:, 4]

#Grafico 3D
fig = plt.figure(1)
ax = fig.gca(projection='3d')
plt.grid()
#titolo e nome assi
ax.set_title('PARTICELLA CARICA IN \n CAMPO ELETTROMAGNETICO', fontsize=20)
ax.set_xlabel('X(t)')
ax.set_ylabel('Y(t)')
ax.set_zlabel('Z(t)')
#limitiamo gli assi altrimenti la loro grandezza non verrebbe delle dimensioni giuste rispetto al plot
ax.set_xlim(np.min(x),np.max(x))
ax.set_ylim(np.min(y),np.max(y))
ax.set_zlim(np.min(z),np.max(z))


#definiamo le variabili per il grafico animato
line, = plt.plot([],[],[],'b-')
dot, = ax.plot([], [], [], 'bo')

#questa è la funzione che anima dove per ogni i metti il valore della soluzione dentro le liste delle variabili sopra definite
def animate(i):
    dot.set_data_3d(x[i],y[i], z[i])
    line.set_data_3d(x[:i],y[:i], z[:i])
    return  dot, line


anim = animation.FuncAnimation(fig, animate, frames=T, interval=1, blit=True, repeat=True) #questa è l'animazione vera e propria, a cui passi il grafico, la funzione che anima, la funzione di init, il numero di fotogrammi, bilt lo lasci true perchè l'animazione è migliore, e repeat si spiega da solo

#anim.save('P_EM.mp4', fps=50, extra_args=['-vcodec', 'libx264']) #salva sulla cartella corrente l'animazione
plt.show()
##
fig1 = plt.figure(2)
ax1 = fig1.gca(projection='3d')
plt.grid()

ax1.set_title('PARTICELLA CARICA IN \n CAMPO ELETTROMAGNETICO', fontsize=20)
ax1.set_xlabel('X(t)')
ax1.set_ylabel('Y(t)')
ax1.set_zlabel('Z(t)')
plt.grid()

ax1.plot(x,y,z)

plt.show()
