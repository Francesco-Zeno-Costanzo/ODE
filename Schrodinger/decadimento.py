import time
import numpy as np
import scipy.special as ssp
from scipy.sparse import diags
import  matplotlib.pyplot  as  plt
start_time=time.time()

n=3000
xr=40
xl=0
L=xr-xl
h=(xr-xl)/(n)
tt=np.linspace(0, n, n)
xp=xl+h*tt
rn=10
V0=-2

def f(x):
    x=np.array([])
    for i in np.arange(xl, xr+xr/n, xr/n):
        if i>rn:
            x=np.insert(x, len(x), 150/i)
        elif i>0 and i<rn:
            x=np.insert(x, len(x), V0)
        elif i==0:
            x=np.insert(x, len(x), 100)
    return x

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
m=20
k=5
plt.figure(1)
plt.ylim(V0-0.1, 17)
plt.title("$\psi(x)$ Barriera coulombiana n=%d, %d" %(m,k), fontsize=20)
plt.ylabel('$\psi(x)$', fontsize=15)
plt.xlabel('x', fontsize=15)
plt.annotate('$r_N$',(rn, V0))
x=np.linspace(xl, xr, n)
plt.grid()

plt.errorbar(xp, psi[:,m]+avals[m], fmt='.', markersize=3, label='$\psi(x)_{%d}$ numerica' %m)
plt.errorbar(xp, psi[:,k]+avals[k], fmt='.', markersize=3, label='$\psi(x)_{%d}$ numerica' %k)
plt.plot(x,f(x), color='red', label='V(x)')
plt.plot(x,np.ones(len(x))*avals[m], color='black', linestyle='--', label='$E_{%d}=%f$' %(m, avals[m]))
plt.plot(x,np.ones(len(x))*avals[k], color='black', linestyle='--', label='$E_{%d}=%f$' %(k, avals[k]))

plt.legend(loc='best')
plt.show()
print("--- %s seconds ---" % (time.time() - start_time-b))
print(avals[0:m])

