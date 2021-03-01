import time
import numpy as np
from scipy.sparse import diags
import  matplotlib.pyplot  as  plt

start_time=time.time()

n=2000
xr=10
xl=-10
L=xr-xl
h=(xr-xl)/(n)
tt=np.linspace(0, n, n)
xp=xl+h*tt
g=0.05
def U(x):
    return -1/2 * x**2 + (g/2)*x**4 + 1/(8*g)

P=diags([1, -2, 1], [-1, 0, 1], shape=(n, n)).toarray()
V=diags(U(xp), 0, shape=(n, n)).toarray()
H=-(1/(2*h**2))*P+V

a=(time.time() - start_time)
print("--- %s seconds ---" %a)

aval, avec=np.linalg.eig(H)

b=(time.time() - start_time - a)
print("--- %s seconds ---" % b)

avals=np.sort(aval)
avecs=avec[:,aval.argsort()]
psi=avecs/np.sqrt(h)
##
m=0
plt.figure(1)
plt.title("$\psi(x)$ doppia buca (ammoniaca) n=%d" %m, fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)$', fontsize=15)
plt.grid()
plt.ylim(np.min(psi[:,m]+avals[m])-0.2, np.max(psi[:,m]+avals[m])+0.2)

x=np.linspace(xl, xr, n)

plt.errorbar(xp, psi[:,m]+avals[m], fmt='.', label='$\psi(x)$ numerica')
plt.plot(x,U(x), color='red', label='V(x)= $- \dfrac{1}{2} x^2 +  \dfrac{g}{2} x^4 + \dfrac{1}{8g}$')
plt.plot(x,np.ones(len(x))*avals[m], color='black', linestyle='--', label='$E_{%d}=%f$' %(m, avals[m]))
plt.legend(loc='best')


plt.show()
print("--- %s seconds ---" % (time.time() - start_time-b-a))
print(avals[0:m+1])
