import time
import numpy as np
import scipy.special as ssp
from scipy.sparse import diags
import  matplotlib.pyplot  as  plt

start_time=time.time()

def fatt(n):
    if(n == 0) or (n == 1):
        return 1
    else:
        return n*fatt(n-1)

def G(x, m, l):
    return (2/(m**2))*np.sqrt((fatt(m-l-1))/ssp.gamma(m+l))*((2*x/(m))**l)*np.exp(-x/m)*ssp.eval_genlaguerre(m-l-1, 2*l+1, 2*x/m)

n=2000
xr=25
xl=0
L=xr-xl
h=(xr-xl)/(n)
tt=np.linspace(0, n, n+1)
tt[0]=tt[0]+1
xp=xl+h*tt

def f(x):
    return -1/x
P=diags([1, -2, 1], [-1, 0, 1], shape=(n+1, n+1)).toarray()
V=diags(f(xp), 0, shape=(n+1, n+1)).toarray()
H=-(1/(2*h**2))*P+V

a= (time.time() - start_time)
print("--- %s seconds ---" %a)

aval, avec=np.linalg.eig(H)

b= (time.time() - start_time - a)
print("--- %s seconds ---" %b)

avals=np.sort(aval)
avecs=avec[:,aval.argsort()]
psi=avecs/np.sqrt(h)
##
m=1
plt.figure(1)
plt.title("$\psi(x)$ Radiale atomo di igrogeno", fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)$', fontsize=15)
plt.errorbar(xp, (psi[:,m-1])**2+avals[m-1], fmt='.',  label='$\psi(x)$ numerica %ds' %m)
plt.plot(x,np.ones(len(x))*avals[m-1], color='black', linestyle='--', label='$E_{%d}=%f$' %(m-1, avals[m-1]))

x=np.linspace(xl+1e-12, xr, n)
plt.plot(x ,f(x), color='black', label='V(x)=$\dfrac{1}{x}$')
plt.ylim(-1,0.75)

x=np.linspace(xl, xr, n+1)
plt.plot(x , (x*G(x, m, 0))**2+avals[m-1], color='red', label='$\psi(x)$ analitica %ds' %m)
plt.grid()
plt.legend(loc='best')

plt.figure(2)
plt.title("Differenza fra $\psi(x)$ numerica e esatta", fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)_{num}-\psi(x)_{es}$', fontsize=15)
plt.errorbar(xp, psi[:,m-1]**2-(x*G(x, m, 0))**2, fmt='.')
plt.grid()

plt.show()
print("--- %s seconds ---" % (time.time() - start_time - b-a))
print(avals[0:m])
