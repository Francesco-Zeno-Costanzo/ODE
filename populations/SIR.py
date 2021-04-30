import numpy as np
import matplotlib.pyplot as plt


N=100
def f(x, y, m=1):
    x_dot=vx
    y_dot=vy
    vx_dot= m * x/np.sqrt(x**2+y**2)
    vy_dot= m * y/np.sqrt(x**2+y**2)
    return x_dot, y_dot, vx_dot, vy_dot


dt = 0.001
num_steps = 10000


xs = np.zeros(num_steps + 1)
ys = np.zeros(num_steps + 1)
vx = np.zeros(num_steps + 1)
vy = np.zeros(num_steps + 1)

xs[0], ys[0], vx[0], vy[0] = (1, 1, 0, 1)


for i in range(num_steps):
    Sk1, Ik1, Rk1 = f(Ss[i], Is[i], Rs[i])
    Sk2, Ik2, Rk2 = f(Ss[i] + Sk1*dt/2, Is[i]+ Ik1*dt/2, Rs[i]+ Rk1*dt/2)
    Sk3, Ik3, Rk3 = f(Ss[i] + Sk2*dt/2, Is[i]+ Ik2*dt/2, Rs[i]+ Rk2*dt/2)
    Sk4, Ik4, Rk4 = f(Ss[i] + Sk3*dt, Is[i]+ Ik3*dt, Rs[i]+ Rk3*dt)
    Ss[i + 1] = Ss[i] + (dt/6)*(Sk1 + 2*Sk2 + 2*Sk3 + Sk4)
    Is[i + 1] = Is[i] + (dt/6)*(Ik1 + 2*Ik2 + 2*Ik3 + Ik4)
    Rs[i + 1] = Rs[i] + (dt/6)*(Rk1 + 2*Rk2 + 2*Rk3 + Rk4)
    ts[i + 1] = ts[i] + dt


plt.figure(1)
plt.title('Modello SIR', fontsize=20)
plt.plot(ts, Ss, 'blue' , label='Suscettibili')
plt.plot(ts, Is, 'red'  , label='Infetti')
plt.plot(ts, Rs, 'black', label='Rimessi')
plt.legend(loc='best')
plt.grid()

plt.show()