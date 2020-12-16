import numpy as np
import matplotlib.pyplot as plt


N=100
def f(S, I, R, c=0.04, g=0.4):
    S_dot = - S*I*c
    I_dot = c*S*I -g*I
    R_dot = g*I
    return S_dot, I_dot, R_dot


dt = 0.001
num_steps = 10000


Ss = np.zeros(num_steps + 1)
Is = np.zeros(num_steps + 1)
Rs = np.zeros(num_steps + 1)
ts = np.zeros(num_steps + 1)

Ss[0], Is[0], Rs[0]=(N-1, 1, 0)


for i in range(num_steps):
    Sk1, Ik1, Rk1 = f(Ss[i], Is[i], Rs[i])
    Sk2, Ik2, Rk2 = f(Ss[i] + Sk1*dt/2, Is[i]+ Ik1*dt/2, Rs[i]+ Rk1*dt/2)
    Sk3, Ik3, Rk3 = f(Ss[i] + Sk2*dt/2, Is[i]+ Ik2*dt/2, Rs[i]+ Rk2*dt/2)
    Sk4, Ik4, Rk4 = f(Ss[i] + Sk3*dt, Is[i]+ Ik3*dt, Rs[i]+ Rk3*dt)
    Ss[i + 1] = Ss[i] + (dt/6)*(Sk1 + 2*Sk2 + 2*Sk3 + Sk4)
    Is[i + 1] = Is[i] + (dt/6)*(Ik1 + 2*Ik2 + 2*Ik3 + Ik4)
    Rs[i + 1] = Rs[i] + (dt/6)*(Rk1 + 2*Rk2 + 2*Rk3 + Rk4)
    ts[i + 1 ]= ts[i] + dt


plt.figure(1)
plt.title('Modello SIR')
plt.plot(ts, Ss, 'blue' , label='Suscettibili')
plt.plot(ts, Is, 'red'  , label='Infetti')
plt.plot(ts, Rs, 'black', label='Rimessi')
plt.legend(loc='best')
plt.grid()

plt.show()