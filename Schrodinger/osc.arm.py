import time
import numpy as np
import scipy.special as ssp
from scipy.sparse import diags
import  matplotlib.pyplot  as  plt

start_time=time.time()
def G(x, m):
    return (1/(np.pi)**(1/4))*(1/np.sqrt((2**m)*ssp.gamma(m+1)))*ssp.eval_hermite(m, x)*np.exp(-(x**2)/2)

n=1000
xr=10
xl=-10
L=xr-xl
h=(xr-xl)/(n)
tt=np.linspace(0, n, n)
xp=xl+h*tt

def U(x):
    return 1/2 * x**2
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
m=10
plt.figure(1)
plt.title("$\psi(x)$ dell'oscillatore armonico n=%d" %m, fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)$', fontsize=15)
plt.grid()
plt.ylim(m-0.2,m+1.3)
x=np.linspace(xl, xr, n)

plt.errorbar(xp, psi[:,m]+avals[m], fmt='.', label='$\psi(x)$ numerica')
plt.plot(x,G(x, m)+avals[m], color='black', label='$\psi(x)$ analitica')
plt.plot(x,U(x), color='red', label='V(x)= $ \dfrac{1}{2} x^2 $')
plt.plot(x,np.ones(len(x))*avals[m], color='black', linestyle='--', label='$E_{%d}=%f$' %(m, avals[m]))
plt.legend(loc='best')

plt.figure(3)
plt.title("Differenza fra $\psi(x)$ numerica e esatta", fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)_{num}-\psi(x)_{es}$', fontsize=15)
plt.grid()

plt.errorbar(xp, psi[:,m]-G(x, m), fmt='.')

plt.show()
print("--- %s seconds ---" % (time.time() - start_time-b-a))
print(avals[0:m]-(1/2 + np.linspace(0, m-1, m)))
