import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt

xi0=10
nP=1000
xiMax=30
step=0.005
DeltaXi=xi0/nP
tau=1.0e-12
xi=np.linspace(0, xiMax, nP)
x=np.linspace(-xiMax,xiMax,(2*nP)+1)
n=3 #stato
EigvStart=-0.7 #circa energia che ci aspettiamo

V0=-1
params=[EigvStart , xi0]
psi=[]

def V():
    x=np.array([])
    for i in np.arange(-xiMax, xiMax + xiMax/nP, xiMax/nP):
        if i>xi0:
            x=np.insert(x, len(x), 0)
        elif i>-xi0 and i<xi0:
            x=np.insert(x, len(x), V0)
        elif i<xi0:
            x=np.insert(x, len(x), 0)
    return x

def F(y, xi, params) :
    psi , dpsidxi=y
    E, xi0=params
    if xi>xi0:
        dydt=[dpsidxi, -E*psi ]
    else :
        dydt=[dpsidxi, (V0-E)*psi ]
    return dydt

def Epsi(params, xi, n, EigvStart, step, tau, F, psi):
    #primo estremo
    E1=EigvStart
    params[0]=E1
    if n%2==0:
        y=[1.0, 0.0]
    else:
        y=[0.0, 1.0]
    sol=si.odeint(F, y, xi, args=(params, ))
    PsiEnd1=sol[-1, 0]
    #secondo estremo
    while True:
        E2=E1+step
        if E2>0: #fuori dalla buca
            return -1
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
        if abs(E1-E2)<tau :
            break
        params[0]=Em
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


eigv=Epsi(params, xi, n ,EigvStart, step, tau, F, psi)
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
psi=np.square(psi) #modulo quadro
psi=[x+eigv for x in psi ] # traslo per comodità grafica

plt.figure(1)
plt.title("$\psi(x)^2$ Buca quadrata n=%d" %n, fontsize=20)
plt.plot(x, np.ones(len(x))*eigv, color='black', linestyle='--', label='$E_{%d}=%f$' %(n, eigv))
plt.plot(x,  psi,label='$\psi(x)$ numerica' )
plt.ylim(V0-0.1, 0.2)
plt.grid()
plt.plot(x, V(), color='black', label='V(x)')
plt.annotate('$V_0$',(xi0, V0))
plt.legend(loc='best')
plt.ylabel('$\psi(x)^2$', fontsize=15)
plt.xlabel('x', fontsize=15)
plt.show()
