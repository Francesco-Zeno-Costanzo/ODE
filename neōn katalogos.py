import numpy as np
import scipy.integrate
import  matplotlib.pyplot  as  plt
o0 = 9
v0 = 0
x0 = 1
xf = 10
'''
x''(t)=-o0*x(t)
x'(0)=v0
x(0)=x0
'''

##Soluzione analitica
def Sol(t):
    return v0/np.sqrt(o0) * np.sin(np.sqrt(o0)*t) + x0*np.cos(np.sqrt(o0)*t)

##Soluzione numerica di odeint
def osc(y, t):
    theta, omega = y
    dydt = [omega,  - o0*theta]
    return dydt

y0 = [x0 , v0] #x(0), x(0)'
t = np.linspace(0, xf, 10001)
sol = scipy.integrate.odeint(osc, y0, t)

##Soluzione numerica con il metodo di eluero
def h(x, v):
    x_dot = v
    v_dot = - o0*x
    return x_dot, v_dot

num_steps = 10000
dt = xf/num_steps


xs = np.zeros(num_steps + 1)
vs = np.zeros(num_steps + 1)
ts = np.zeros(num_steps + 1)

xs[0], vs[0]=(x0, v0)

for i in range(num_steps):
    x1, v1 = h(xs[i], vs[i])
    xs[i + 1] = xs[i] + dt*x1
    vs[i + 1] = vs[i] + dt*v1
    ts[i + 1] = ts[i] + dt

##Soluzione numerica con il metodo di eluero semi implicito (integratore simplettico)
def f(v):
    x_dot = v
    return x_dot
def g(x):
    v_dot = - o0*x
    return v_dot

num_steps = 10000
dt = xf/num_steps

xs1 = np.zeros(num_steps + 1)
vs1 = np.zeros(num_steps + 1)
ts1 = np.zeros(num_steps + 1)

xs1[0], vs1[0]=(x0, v0)

for i in range(num_steps):
    vs1[i + 1] = vs1[i] + dt*g(xs1[i])
    xs1[i + 1] = xs1[i] + dt*f(vs1[i+1])
    ts1[i + 1] = ts1[i] + dt

##Soluzione numerica con il metodo di eluero implicito
'''
y'=f(t, y)
y[i+1]=y[i]+dt*f(t[i+1], y[i+1])
nel nostro caso va risolto il sitema
v[i+1] = v[i] -o0*dt*x[i+1]
x[i+1] = x[i]    +dt*v[i+1]
Se il sistema non è risolubile analiticamente si può usare un metodo numerico.
Un esempio è al metodo del punto medio implicito
'''
num_steps = 10000
dt = xf/num_steps

xs2 = np.zeros(num_steps + 1)
vs2 = np.zeros(num_steps + 1)
ts2 = np.zeros(num_steps + 1)

xs2[0], vs2[0]=(x0, v0)

for i in range(num_steps):
    vs2[i + 1] = (vs2[i] -o0*dt*xs2[i])/(1+o0*dt**2)
    xs2[i + 1] = xs2[i] + dt*vs2[i+1]
    ts2[i + 1] = ts2[i] + dt

##Soluzione numerica con il metodo velocity verlet (integratore simplettico)
def F(x): #forza che agisce sul sistema
    return - o0*x

num_steps = 10000
dt = xf/num_steps

xs3 = np.zeros(num_steps+1)
vs3 = np.zeros(num_steps+1)
ts3 = np.zeros(num_steps+1)

xs3[0], vs3[0] = (x0, v0)

for i in range(num_steps):
    Fx = F(xs3[i])
    xs3[i + 1] = xs3[i] + vs3[i]*dt + Fx*((dt**2)/2)
    F1x= F(xs3[i+1])
    vs3[i + 1] = vs3[i] + (Fx+F1x)*dt/2
    ts3[i + 1] = ts3[i] + dt

##Soluzione numerica con il metodo runge-kutta di ordine 4
def r(x, v):
    x_dot = v
    v_dot = - o0*x
    return x_dot, v_dot

num_steps = 10000
dt = xf/num_steps

xs4 = np.zeros(num_steps + 1)
vs4 = np.zeros(num_steps + 1)
ts4 = np.zeros(num_steps + 1)

xs4[0], vs4[0]=(x0, v0)


for i in range(num_steps):
    xk1, vk1 = r(xs4[i], vs4[i])
    xk2, vk2 = r(xs4[i] + xk1*dt/2, vs4[i]+ vk1*dt/2)
    xk3, vk3 = r(xs4[i] + xk2*dt/2, vs4[i]+ vk2*dt/2)
    xk4, vk4 = r(xs4[i] + xk3*dt, vs4[i]+ vk3*dt)
    xs4[i + 1] = xs4[i] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
    vs4[i + 1] = vs4[i] + (dt/6)*(vk1 + 2*vk2 + 2*vk3 + vk4)
    ts4[i + 1] = ts4[i] + dt

##Soluzione numerica con il metodo del punto medio implicito
'''
y'=f(t, y)
y[i+1]=y[i]+dt*f((t[i+1]+t[i])/2, (y[i+1]+y[i])/2)
nel nostro caso va risolto il sitema
v[i+1] = v[i] -o0*dt*(x[i+1] + x[i])/2
x[i+1] = x[i]    +dt*(v[i+1] + v[i])/2
Se in sistema non è risolubile analiticamente possiamo procedere numericamente.
Il seguente è un esempio per x''=-sin(o0*x):

import scipy.optimize as so
def F(V, x0, v0):
    v1, x1= V
    R1=v1-v0-dt*(-o0*np.sin((x1+x0)/2))
    R2=x1-x0-dt*(v1+v0)/2
    return [R1, R2]

for i in range(num_steps):
    xstart=(ys[i], xs[i])
    ys[i + 1], xs[i + 1] = so.fsolve(F , xstart, args=(xs[i], ys[i]))
    ts[i + 1] = ts[i] + dt
#dove xstart dipende da i per velocizzare, su grandi tempi, il calcolo
'''

num_steps = 10000
dt = xf/num_steps

xs5 = np.zeros(num_steps + 1)
vs5 = np.zeros(num_steps + 1)
ts5 = np.zeros(num_steps + 1)

xs5[0], vs5[0]=(x0, v0)

for i in range(num_steps):
    vs5[i + 1] = (vs5[i] -o0*dt*xs5[i] - o0*(dt**2)*vs5[i]/4)/(1+(o0*dt**2)/4)
    xs5[i + 1] = xs5[i] + dt*vs5[i]/2 + dt*vs5[i+1]/2
    ts5[i + 1] = ts5[i] + dt

##Soluzione numerica con il metodo ruth3-4 (integratore simplettico)

def F(x):
    return - o0*x

#coefficienti ruth 3
c = np.array([2/3, -2/3, 1])
d = np.array([7/24, 3/4, -1/24])
'''
#coefficienti ruth 4
l=2**(1/3)
c=np.array([1/(2*(2-l)), (1-l)/(2*(2-l)), (1-l)/(2*(2-l)), 1/(2*(2-l))])
d=np.array([0, 1/(2-l), -l/(2-l), 1/(2-l)])
'''

num_steps = 10000
dt = xf/num_steps

xs6 = np.zeros(num_steps+1)
vs6 = np.zeros(num_steps+1)
ts6 = np.zeros(num_steps+1)

xs6[0], vs6[0] = (x0, v0)

for i in range(num_steps):
    x, v = xs6[i], vs6[i]
    for k in range(len(d)):
        v+=d[k]*F(x)*dt
        x+=c[k]*v*dt
    xs6[i + 1] = x
    vs6[i + 1] = v
    ts6[i + 1] = ts6[i] + dt

##Grafico soluzioni
plt.figure(1)
plt.title('Confronto soluzioni', fontsize=20)
plt.xlabel('t', fontsize=15)
plt.ylabel(r'$\vartheta(t)$', fontsize=15)
plt.plot(t, Sol(t), 'black', label='sol analitica')
plt.plot(t, sol[:, 0], 'blue', label='odeint')
plt.plot(ts, xs, 'red', label='Eulero')
plt.plot(ts1, xs1, 'green', label='Eulero semi implicito (integratore simplettico)')
plt.plot(ts2, xs2, 'brown', label='Eulero implicito I')
plt.plot(ts3, xs3, 'yellow', label='velocity verlet (integratore simplettico)')
plt.plot(ts4, xs4, 'pink', label='Runge Kutta 4')
plt.plot(ts5, xs5, 'orange', label='punto medio implicito (integratore simplettico)')
plt.plot(ts6, xs6, 'fuchsia', label='ruth (integratore simplettico)')
plt.legend(loc='best')
plt.grid()

##Grafico spazio delle fasi
plt.figure(2)
plt.suptitle('Spazio delle fasi', fontsize=20)
plt.grid()
plt.subplot(421)
plt.plot(sol[:, 0], sol[:, 1], 'k', label='odeint')
plt.legend(loc='best')
plt.grid()
plt.subplot(422)
plt.plot(xs, vs, 'k', label='Eulero')
plt.legend(loc='best')
plt.grid()
plt.subplot(423)
plt.plot(xs1, vs1, 'k', label='Eulero semi implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(424)
plt.plot(xs2, vs2, 'k', label='Eulero implicito')
plt.legend(loc='best')
plt.grid()
plt.subplot(425)
plt.plot(xs3, vs3, 'k', label='velocity verlet (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(426)
plt.plot(xs4, vs4, 'k', label='Runge Kutta 4')
plt.legend(loc='best')
plt.grid()
plt.subplot(427)
plt.plot(xs5, vs5, 'k', label='punto medio implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(428)
plt.plot(xs6, vs6, 'k', label='ruth (integratore simplettico)')
plt.legend(loc='best')
plt.grid()

##Grafico differenze
plt.figure(3)
plt.subplot(421)
plt.title('Differenza tra soluzione esatta e numerica', fontsize=20)
plt.plot(t, Sol(t)-sol[:, 0], 'k', label='odeint')
plt.legend(loc='best')
plt.grid()
plt.subplot(422)
plt.plot(t, Sol(t)-xs, 'k', label='Eulero')
plt.legend(loc='best')
plt.grid()
plt.subplot(423)
plt.plot(t, Sol(t)-xs1, 'k', label='Eulero semi implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(424)
plt.plot(t, Sol(t)-xs2, 'k', label='Eulero implicito')
plt.legend(loc='best')
plt.grid()
plt.subplot(425)
plt.plot(t, Sol(t)-xs3, 'k', label='velocity verlet (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(426)
plt.plot(t, Sol(t)-xs4, 'k', label='Runge Kutta 4')
plt.legend(loc='best')
plt.grid()
plt.subplot(427)
plt.plot(t, Sol(t)-xs5, 'k', label='punto medio implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(428)
plt.plot(t, Sol(t)-xs6, 'k', label='ruth (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
##Grafico dell'energia
def U(v, x):
    return (v**2 +o0*x**2)-(v[0]**2 +o0*x[0]**2)

plt.figure(4)
plt.subplot(421)
plt.title('Differenza fra enerigia iniziale', fontsize=20)
plt.plot(t, U(sol[:, 1], sol[:,0]) , 'k', label='odeint')
plt.legend(loc='best')
plt.grid()
plt.subplot(422)
plt.title(' ed energia al tempo t del sistema', fontsize=20)
plt.plot(t, U(vs, xs), 'k', label='Eulero')
plt.legend(loc='best')
plt.grid()
plt.subplot(423)
plt.plot(t, U(vs1, xs1), 'k', label='Eulero semi implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(424)
plt.plot(t, U(vs2, xs2), 'k', label='Eulero implicito')
plt.legend(loc='best')
plt.grid()
plt.subplot(425)
plt.plot(t, U(vs3, xs3), 'k', label='velocity verlet (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(426)
plt.plot(t, U(vs4, xs4), 'k', label='Runge Kutta 4')
plt.legend(loc='best')
plt.grid()
plt.subplot(427)
plt.plot(t, U(vs5, xs5), 'k', label='punto medio implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(428)
plt.plot(t, U(vs6, xs6), 'k', label='ruth (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.show()
