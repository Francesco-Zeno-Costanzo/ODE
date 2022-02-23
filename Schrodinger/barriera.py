import time
import numpy as np
import scipy.special as ssp
from scipy.sparse import diags
import  matplotlib.pyplot  as  plt

start_time=time.time()

n=2000
xr=15
xl=-15
L=xr-xl
h=(xr-xl)/(n)
tt=np.linspace(0, n, n)
xp=xl+h*tt

A=50
s=2
def f(x):
    return A*np.exp(-((x/s)**2))
P=diags([1, -2, 1], [-1, 0, 1], shape=(n, n)).toarray()
V=diags(f(xp), 0, shape=(n, n)).toarray()
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
m=40
plt.figure(1)
plt.ylim(np.min(psi[:,m]+avals[m])-1/2, np.max(psi[:,m]+avals[m])+1/2)
plt.title("$\psi(x)$ Barriera n=%d" %m, fontsize=20)
plt.ylabel('$\psi(x)$', fontsize=15)
plt.xlabel('x', fontsize=15)
x=np.linspace(xl, xr, n)
plt.grid()

plt.errorbar(xp, psi[:,m]+avals[m], fmt='.', label='$\psi(x)$ numerica')
plt.plot(x,f(x), color='red', label='V(x)=$Ae^{-(x / \sigma)^2}$')
plt.plot(x,np.ones(len(x))*avals[m], color='black', linestyle='--', label='$E_{%d}=%f$' %(m, avals[m]))

plt.legend(loc='best')
plt.show()

c=(time.time() - start_time-b-a)
print("--- %s seconds ---" % c)


plt.show()
