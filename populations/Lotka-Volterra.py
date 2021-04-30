import numpy as np
import scipy.special
import scipy.integrate
import  matplotlib.pyplot  as  plt


def f(y, t, a11, a12, a21, a22, b1, b2, c1, c2):
    x, z = y
    dydt = [b1*x + c1*z + (a11-a12*z)*x, c2*x + b2*z + (a21*x-a22)*z]
    return dydt

a11=1
a12=0.2
a21=0.2
a22=0.8
b1=0
b2=0
c1=0
c2=0
x0 = [10, 1]     # x(0), x(0)',
t = np.linspace(0, 20, 1000)
sol = scipy.integrate.odeint(f, x0, t, args=( a11, a12, a21, a22, b1, b2, c1, c2))

plt.figure(1)
plt.title('Sistema Lotka-Volterra')
plt.xlabel('t')
plt.ylabel('popolazione')

plt.plot(t, sol[:,0], color='green', label='prede')
plt.plot(t, sol[:,1], color='red', label='predatori')

plt.grid()
plt.legend(loc='best')

plt.figure(2)
plt.title('Spazio delle fasi')
plt.plot(sol[:,1], sol[:,0])
plt.grid()


plt.show()
