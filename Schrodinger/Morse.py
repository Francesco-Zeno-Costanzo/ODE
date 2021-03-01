import time
import numpy as np
import scipy.special as ssp
from scipy.sparse import diags
import  matplotlib.pyplot  as  plt
start_time=time.time()

n=3000
xr=10
xl=-0
L=xr-xl
h=(xr-xl)/(n)
tt=np.linspace(0, n, n)
xp=xl+h*tt

x0=1
B=152
def f(x):
    return B*(1 - np.exp(-(x-x0)))**2
l=np.sqrt(2*B)
def G(x, n):
    return np.sqrt(((ssp.gamma(n+1)*(2*l-2*n-1))/ssp.gamma(2*l-n)))*(2*l*np.exp(-(x-x0)))**(l-n-1/2)*np.exp(-1/2 *2*l*np.exp(-(x-x0)))*ssp.eval_genlaguerre(n, 2*l-2*n-1, 2*l*np.exp(-(x-x0)))

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
m=14
x=np.linspace(xl, xr, n)
plt.figure(1)
plt.ylim(np.min(psi[:,m]+avals[m]), np.max(psi[:,m]+avals[m]))
plt.title("$\psi(x)$ Potenziale di Morse n=%d" %m, fontsize=20)
plt.ylabel('$\psi(x)$', fontsize=15)
plt.xlabel('x', fontsize=15)
plt.grid()

plt.errorbar(xp, psi[:,m]+avals[m], fmt='.', label='$\psi(x)$ numerica')
plt.plot(x, G(x, m)+avals[m], color='red', label='$\psi(x)$ analitica')
plt.plot(x, f(x), color='black', label='V(x)= $B(1-e^{-(x-x_0)})^2$')
plt.plot(x,np.ones(len(x))*avals[m], color='black', linestyle='--', label='$E_{%d}=%f$' %(m, avals[m]))
plt.legend(loc='best')


plt.figure(3)
plt.title("Differenza fra $\psi(x)$ numerica e esatta", fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)_{num}-\psi(x)_{es}$', fontsize=15)
plt.errorbar(xp, psi[:,m]-G(x, m), fmt='.')
plt.grid()

plt.show()
print("--- %s seconds ---" % (time.time() - start_time-b-a))

E=np.array([])
N=int(np.sqrt(2*B)-0.5)
for i in range(N):
    o=np.sqrt(2*B)
    En=o*(i+1/2)-((o**2)/(4*B))*(i+1/2)**2
    E=np.insert(E, len(E), En)
print(avals[0:N]-E)


