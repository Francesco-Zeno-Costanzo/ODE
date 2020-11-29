import numpy as np
import scipy.special
import scipy.integrate
import  matplotlib.pyplot  as  plt


def pend(y, t, b, c):
    theta, omega = y
    dydt = [omega, -b*omega  - c*np.sin(theta)]
    return dydt

b=0.25
c=5

y0 = [np.pi-3, 0] #x(0), x(0)'
t = np.linspace(0, 60, 1000)
sol = scipy.integrate.odeint(pend, y0, t, args=(b, c))

def osc(l, t, b, c):
    theta, omega = l
    dldt = [omega, -b*omega  - c*theta]
    return dldt

b=0.25
c=5

y1 = [np.pi-3 , 0] #x(0), x(0)'
sol1 = scipy.integrate.odeint(osc, y1, t, args=(b, c))

plt.figure(1)
plt.title('Confronto fra pendolo e oscillatore', fontsize=20)
plt.xlabel('t')
plt.ylabel('$\phi(t)$')

plt.plot(t, sol[:, 0], 'b', label='$\phi(t)$ pendolo')
plt.plot(t, sol1[:, 0], 'g', label='$\phi(t)$ oscillatore')

plt.legend(loc='best')
plt.grid()

plt.figure(2)
plt.title('$\phi(t)_{pendolo} - \phi(t)_{oscillatore}$', fontsize=20)
plt.xlabel('t')
plt.ylabel('$\phi(t)_{pendolo} - \phi(t)_{oscillatore}$')

plt.plot(t, sol[:,0]-sol1[:,0])

plt.grid()
plt.show()