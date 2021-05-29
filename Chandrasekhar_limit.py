import time
import numpy as np
import matplotlib.pyplot as plt
from scipy import constants as cst

#importiamo costanti varie ed eventuali, ci si fida di scipy e di wiki
c=cst.value(u'speed of light in vacuum')
h=cst.value(u'Planck constant')
G=cst.value(u'Newtonian constant of gravitation')
mp=cst.value(u'proton mass')
MS=1.98892*1e30
mi=2

#equazione da risolvere
def LE(phi, psi, x, n):

    phi_dot = psi
    psi_dot = -phi**n - 2*psi/x

    return phi_dot, psi_dot

#integratore
def RK4(n):
    num_steps = 2000000
    dx = 8/num_steps

    phi = np.zeros(num_steps + 1)
    psi = np.zeros(num_steps + 1)
    x = np.zeros(num_steps + 1)

    phi[0], psi[0]=(1, 0) #condizioni in modo che la densità al centro sia effettivamente lei
    x[0]=1e-10 #evitiamo divergenze

    for i in range(num_steps):

        pk1, sk1 = LE(phi[i], psi[i], x[i], n)
        pk2, sk2 = LE(phi[i] + pk1*dx/2, psi[i] + sk1*dx/2, x[i] + dx/2, n)
        pk3, sk3 = LE(phi[i] + pk2*dx/2, psi[i] + sk2*dx/2, x[i] + dx/2, n)
        pk4, sk4 = LE(phi[i] + pk3*dx, psi[i] + sk3*dx, x[i] + dx, n)

        phi[i + 1] = phi[i] + (dx/6)*(pk1 + 2*pk2 + 2*pk3 + pk4)
        psi[i + 1] = psi[i] + (dx/6)*(sk1 + 2*sk2 + 2*sk3 + sk4)
        x[i + 1] = x[i] + dx

    return phi, psi, x

#soluzioni analitiche
def f0(x):
    return 1 - (x**2)/6

def f1(x):
    return np.sin(x)/x

def f5(x):
    return 1/np.sqrt(1 + (x**2)/3)

#risolviamo le equazioni anche nei casi analitici per vedere se l'integratore si comporta bene
start_time=time.time()

phi0, psi0, w = RK4(0)
phi1, psi1, w = RK4(1)
phi5, psi5, w = RK4(5)
phi3, psi3, x = RK4(3)

A=(time.time() - start_time)
print(" --- %s secondi per la risoluzione --- " %A)

#cerchiamo lo zero con il metodo di bisezione
a=0
b=len(x)-1
fa=phi3[a]
fb=phi3[b]

iter=1

tau = x[1] - x[0] #chiaramente una precisione maggiore creerebbe un loop infinito

while (b-a)>tau:
    d=(a+b)/2.0
    fd=phi3[int(d)]
    if fd==0:
        print("x = " ,d)
    if fd*fa>0:
        a=d
        fa=fd
    else:
        b=d
        fb=fd
    iter+=1

B=(time.time() - start_time - A)
print("--- %s secondi la ricerca dello zero ---" %B)
print('ovvero un numero di interazioni pari a: %d' %iter)

#calcoliamo il necessario per ottenere la massa limite
x0=x[int(d)]
H=psi3[int(d)]*x0**2
MCH=-H*(np.sqrt(6)/(8*np.pi))*(((h*c)/G)**(3/2))*(1/(mp*mi))**2

print('zero della soluzione: %f' %x0)
print("x^2*phi'calcolato in x0: %f" %H)
print('massa limite in unità di masse solari: %f' %(MCH/MS))

l=np.linspace(0.1, np.min( psi3*x**2), 1000)
k=np.linspace(0, x0, 1000)
ll=np.ones(len(l))

plt.figure(1)
plt.subplot(311)
plt.title('Equazione di Lane Emden', fontsize=20)
plt.grid()
plt.ylabel('n=0', fontsize=20)
plt.plot(w, f0(w), label='$\phi_0$ analitica')
plt.plot(w, phi0,  label='$\phi_0$ numerica')
plt.legend(loc='best')

plt.subplot(312)
plt.grid()
plt.ylabel('n=1', fontsize=20)
plt.plot(w, f1(w), label='$\phi_1$ analitica')
plt.plot(w, phi1,  label='$\phi_1$ numerica')
plt.legend(loc='best')

plt.subplot(313)
plt.ylabel('n=5', fontsize=20)
plt.plot(w, f5(w), label='$\phi_5$ analitica')
plt.plot(w, phi5,  label='$\phi_5$ numerica')
plt.legend(loc='best')
plt.xlabel(r'$x= \frac{r}{\lambda}$', fontsize=15)
plt.grid()

plt.figure(2)
plt.subplot(311)
plt.title('Differenze fra soluzioni analitiche e numeriche', fontsize=20)
plt.grid()
plt.ylabel('n=0', fontsize=20)
plt.plot(w, f0(w)-phi0)

plt.subplot(312)
plt.ylabel('n=1', fontsize=20)
plt.grid()
plt.plot(w, f1(w)-phi1)
plt.subplot(313)

plt.grid()
plt.ylabel('n=5', fontsize=20)
plt.plot(w, f5(w)-phi5)
plt.xlabel(r'$x= \frac{r}{\lambda}$', fontsize=15)


plt.figure(3)
plt.title('Equazione di Lane Emden per n=3', fontsize=20)
plt.xlabel(r'$x= \frac{r}{\lambda}$', fontsize=15)
plt.plot(x, phi3, label=r'$\phi(t)$')
plt.plot(x, psi3*x**2, label=r'$x^2 \frac{d \phi}{dx}$')
plt.annotate(str(round(x0, 3)), (x0+0.01, l[-1]+0.01))
plt.annotate(str(round(H, 3)), (0, H+0.01))
plt.plot(x0*ll, l, color='black', linestyle='--')
plt.plot(k, H*ll, color='black', linestyle='--')
plt.legend(loc='best')
plt.grid()
plt.show()
