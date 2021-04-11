import numpy as np
import  matplotlib.pyplot  as  plt
'''
problema con condizioni al bordo
x''(t)=3/2 x(t)
x(0)=x0
x(1)=x1
si consiglia l'esecuzione a celle per avere un'idea del valore iniziale da dare
we recommend the execution in cells to get an idea of ​​the initial value to be given
'''
xi=0        #estremo sinistro dell'intervallo
xf=1        #estremo desto dell'intervallo
N=1000      #numero di punti
x0=4        #valore all'estremo sinistro
x1=1        #valore che vogliamo assumi la soluzione nell'estemo destro
tau=1.0e-10 #tolleranza


def f(x, v):
    x_dot = v
    v_dot = (3/2)*x**2
    return x_dot, v_dot

def Rk4(xf, N, x0, s):
    dt = xf/N
    xs = np.zeros(N + 1)
    vs = np.zeros(N + 1)
    xs[0], vs[0]=(x0, s)
    for i in range(N):
        xk1, vk1 = f(xs[i], vs[i])
        xk2, vk2 = f(xs[i] + xk1*dt/2, vs[i]+ vk1*dt/2)
        xk3, vk3 = f(xs[i] + xk2*dt/2, vs[i]+ vk2*dt/2)
        xk4, vk4 = f(xs[i] + xk3*dt, vs[i]+ vk3*dt)
        xs[i + 1] = xs[i] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
        vs[i + 1] = vs[i] + (dt/6)*(vk1 + 2*vk2 + 2*vk3 + vk4)
    return xs

def F(xf, N, x0, start, step, x1, n):
    P=np.zeros(n)
    a=start
    t=np.linspace(start, start+n*step, n)
    for j in range(n):
        P[j]=Rk4(xf, N, x0, a)[-1]
        a+=step
    return t, P-x1

t, y=F(xf, N, x0, -40, 0.0833, x1, 408)
#per vedere orientativamente dove sono gli zeri
plt.figure(1)
plt.title('$\simeq$ Funzione di cui trovare gli zeri', fontsize=20)
plt.ylabel('x(1;s)-x(1)', fontsize=15)
plt.xlabel('s', fontsize=15)
plt.grid()
plt.plot(t, 0*t, color='red', linestyle='--')
plt.plot(t, y)
plt.show()
##
def SH(xf, N, x0, start, step, x1, tau):
    a=start
    sol=Rk4(xf, N, x0, a)
    k=sol[-1]-x1
    while True:
        b=a+step
        sol=Rk4(xf, N, x0, b)
        D=sol[-1]-x1
        if (k*D)<0.0:
            break
        k=D
        a=b
    while abs(a - b)>tau:
        m=(a + b)/2.0
        sol=Rk4(xf, N, x0, m)
        M=sol[-1]-x1
        if (M*k)>0 :
            k=M
            a=m
        else :
            D=M
            b=m
    return m, sol

t1=np.linspace(xi, xf, N+1)
v0, y1=SH(xf, N, x0, -40, 0.1, x1, tau)
v1, y2=SH(xf, N, x0, -10, 0.1, x1, tau)

plt.figure(2)
plt.title('Soluzioni per due diverse condizioni inziali che danno la stessa condizione al contorno \n $x_1(t_f)$-x1=%e, $x_2(t_f)$-x1=%e ' %(y1[-1]-x1, y2[-1]-x1), fontsize=20)

plt.ylabel('x(t)', fontsize=15)
plt.xlabel('t', fontsize=15)
plt.grid()
plt.plot(t1, y1, 'k', label='$x_1(t)$, $\dot{x}(t=0)$=%f' %v0)
plt.plot(t1, y2, 'b', label='$x_2(t)$, $\dot{x}(t=0)$=%f' %v1)
plt.legend(loc='best')
plt.show()
