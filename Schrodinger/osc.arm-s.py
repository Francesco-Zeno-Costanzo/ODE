import numpy as np
import scipy.special as ssp
import scipy.integrate as si
import matplotlib.pyplot as plt

nP=10000
xiMax=10
step=0.005
DeltaXi=xiMax/nP
tau=1.0e-12
xi=np.linspace(0, xiMax, nP)
x=np.linspace(-xiMax,xiMax,(2*nP)+1)

n=16 #stato
EigvStart=n #energia che ci aspettiamo

params=[EigvStart]
psi=[]

def G(x, m):
    return (1/(np.pi)**(1/4))*(1/np.sqrt((2**m)*ssp.gamma(m+1)))*ssp.eval_hermite(m, x)*np.exp(-(x**2)/2)

def f(x):
    return 1/2 * x**2

def F(y, xi, params) :
    psi , dpsi=y
    E,=params
    dydt=[dpsi , ( xi**2 - 2.0*E)*psi ]
    return dydt

def Epsi(params, xi, n, EigvStart, step, tau, F, psi):
    #primo estremo
    E1=EigvStart
    params[0]=E1
    if n%2==0: #parità della spi
        y=[1.0, 0.0]
    else:
        y=[0.0, 1.0]
    sol=si.odeint(F, y, xi, args=(params, ))
    PsiEnd1=sol[-1,0]
    #secondo estremo
    while True:
        E2=E1+step
        params[0]=E2
        sol=si.odeint(F, y, xi, args=(params, ))
        PsiEnd2=sol[-1, 0]
        if (PsiEnd1*PsiEnd2)<0.0:
            break
        PsiEnd1=PsiEnd2
        E1=E2
    #ricerca dell'autovalore
    while True:
        Em=(E1+E2)/2.0
        params[0]=Em
        if abs(E1-E2)<tau :
            break
        sol=si.odeint(F, y, xi, args=(params, ))
        PsiEndM=sol[-1, 0]
        if (PsiEndM*PsiEnd1)>0 :
            PsiEnd1=PsiEndM
            E1=Em
        else :
            PsiEnd2=PsiEndM
            E2=Em
    del psi[ : ]
    for i in range(len(xi)) :
        psi.append(sol[i, 0])
    return Em


eigv=Epsi(params, xi, n, EigvStart ,step, tau , F, psi)
print("E_%d = %f"%(n, eigv))

# elimino divergenze
while True:
    if abs(psi[-2])>abs( psi[-1]):
        break
    psi.pop()
if len(psi)<(nP+1):
    while len(psi)<(nP+1):
        psi.append(0.0)

#normalizzazione
Norm=np.sqrt(2.0*si.simps(np.square(psi), dx=DeltaXi, even='first' ))
normpsi=psi/Norm
psineg=list(reversed(normpsi))
if n%2==1: #  per n dispari le funzioni sono dispare
    for k in range( len(psineg ) ) :
        psineg[k]=-psineg[k]

#attacco le parti
psineg.pop()
normpsi=list(normpsi)
psi=psineg+normpsi
psi=[x+eigv for x in psi ] # traslo per comodità grafica

plt.figure(1)
plt.title("$\psi(x)$ Buca quadrata n=%d" %n, fontsize=20)
plt.errorbar(x,  psi, fmt='.', label='$\psi(x)$ numerica' )
plt.plot(x, np.ones(len(x))*eigv, color='black', linestyle='--', label='$E_{%d}=%f$' %(n, eigv))
plt.plot(x, G(x, n)+eigv, color='black', label='$\psi(x)$ analitica')
plt.plot(x, f(x), color='red', label='V(x)= $ \dfrac{1}{2} x^2 $')
plt.ylim(n-0.2,n+1.3)
plt.grid()
plt.legend(loc='best')
plt.ylabel('$\psi(x)$', fontsize=15)
plt.xlabel('x', fontsize=15)

plt.figure(2)
plt.title("Differenza fra $\psi(x)$ numerica e esatta", fontsize=20)
plt.xlabel('x', fontsize=15)
plt.ylabel('$\psi(x)_{num}-\psi(x)_{es}$', fontsize=15)
plt.grid()
plt.errorbar(x, (psi-G(x, n))-eigv, fmt='.')
plt.show()
