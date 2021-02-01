import numpy as np
import scipy.integrate
import  matplotlib.pyplot  as  plt
o0=9
v0=1
x0=0
'''
x''(t)=-o*x(t)
x'(0)=v0
x(0)=x0
'''
##soluzione analitica
def Sol(t):
    return v0/np.sqrt(o0) * np.sin(np.sqrt(o0)*t) + x0*np.cos(np.sqrt(o0)*t)

##soluzione numerica di odeint
def osc(y, t):
    theta, omega = y
    dydt = [omega,  - o0*theta]
    return dydt

y0 = [x0 , v0] #x(0), x(0)'
t = np.linspace(0, 10, 10001)
sol = scipy.integrate.odeint(osc, y0, t)

##soluzione numerica con il metodo di eluero
def h(x, v):
    x_dot = v
    v_dot = - o0*x
    return x_dot, v_dot

dt = 0.001
num_steps = 10000

xs = np.zeros(num_steps + 1)
vs = np.zeros(num_steps + 1)
ts = np.zeros(num_steps + 1)

xs[0], vs[0]=(x0, v0)

for i in range(num_steps):
    x1, v1 = h(xs[i], vs[i])
    xs[i + 1] = xs[i] + dt*x1
    vs[i + 1] = vs[i] + dt*v1
    ts[i + 1] = ts[i] + dt

##soluzione numerica con il metodo di eluero semi implicito (integratore simplettico)
def f(v):
    x_dot = v
    return x_dot
def g(x):
    v_dot = - o0*x
    return v_dot

dt = 0.001
num_steps = 10000

xs1 = np.zeros(num_steps + 1)
vs1 = np.zeros(num_steps + 1)
ts1 = np.zeros(num_steps + 1)

xs1[0], vs1[0]=(x0, v0)

for i in range(num_steps):
    vs1[i + 1] = vs1[i] + dt*g(xs1[i])
    xs1[i + 1] = xs1[i] + dt*f(vs1[i+1])
    ts1[i + 1] = ts1[i] + dt

##soluzione numerica con il metodo runge-kutta di ordine 4
def r(x, v):
    x_dot = v
    v_dot = - o0*x
    return x_dot, v_dot

dt = 0.001
num_steps = 10000

xs2 = np.zeros(num_steps + 1)
vs2 = np.zeros(num_steps + 1)
ts2 = np.zeros(num_steps + 1)

xs2[0], vs2[0]=(x0, v0)


for i in range(num_steps):
    xk1, vk1 = r(xs2[i], vs2[i])
    xk2, vk2 = r(xs2[i] + xk1*dt/2, vs2[i]+ vk1*dt/2)
    xk3, vk3 = r(xs2[i] + xk2*dt/2, vs2[i]+ vk2*dt/2)
    xk4, vk4 = r(xs2[i] + xk3*dt, vs2[i]+ vk3*dt)
    xs2[i + 1] = xs2[i] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
    vs2[i + 1] = vs2[i] + (dt/6)*(vk1 + 2*vk2 + 2*vk3 + vk4)
    ts2[i + 1] = ts2[i] + dt

##grafico
plt.figure(1)
plt.plot(t, Sol(t), 'black', label='sol analitica')
plt.plot(t, sol[:, 0], 'blue', label='odeint')
plt.plot(ts, xs, 'red', label='Eulero')
plt.plot(ts1, xs1, 'green', label='Eulero semi implicito (integratore simplettico)')
plt.plot(ts2, xs2, 'yellow', label='Runge Kutta 4')
plt.legend(loc='best')
plt.grid()

##grafico differenza
plt.figure(2)
plt.subplot(411)
plt.title('Differenza tra soluzione esatta e numerica', fontsize=20)
plt.plot(t, Sol(t)-sol[:, 0], 'k', label='odeint')
plt.legend(loc='best')
plt.grid()
plt.subplot(412)
plt.plot(t, Sol(t)-xs, 'k', label='Eulero')
plt.legend(loc='best')
plt.grid()
plt.subplot(413)
plt.plot(t, Sol(t)-xs1, 'k', label='Eulero semi implicito (integratore simplettico)')
plt.legend(loc='best')
plt.grid()
plt.subplot(414)
plt.plot(t, Sol(t)-xs2, 'k', label='Runge Kutta 4')
plt.legend(loc='best')
plt.grid()
plt.show()