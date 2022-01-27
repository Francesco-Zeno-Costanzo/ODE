import numpy as np
import matplotlib.pyplot as plt


S0 = 97
I0 = 3
b = 0.4
g = 0.04
def f(S, I):
    S_dot = -(b*I*S)/S0
    I_dot = (b*I*S)/S0 - g*I
    R_dot = g*I
    return S_dot, I_dot, R_dot


dt = 0.001
num_steps = 100000


S = np.zeros(num_steps + 1)
I = np.zeros(num_steps + 1)
R = np.zeros(num_steps + 1)
t = np.zeros(num_steps + 1)

S[0], I[0] = (S0, I0)


for i in range(num_steps):
    Sk1, Ik1, Rk1 = f(S[i], I[i])
    Sk2, Ik2, Rk2 = f(S[i] + Sk1*dt/2, I[i] + Ik1*dt/2)
    Sk3, Ik3, Rk3 = f(S[i] + Sk2*dt/2, I[i] + Ik2*dt/2)
    Sk4, Ik4, Rk4 = f(S[i] + Sk3*dt, I[i] + Ik3*dt)
    S[i + 1] = S[i] + (dt/6)*(Sk1 + 2*Sk2 + 2*Sk3 + Sk4)
    I[i + 1] = I[i] + (dt/6)*(Ik1 + 2*Ik2 + 2*Ik3 + Ik4)
    R[i + 1] = R[i] + (dt/6)*(Rk1 + 2*Rk2 + 2*Rk3 + Rk4)
    t[i + 1] = t[i] + dt


plt.figure(1)
plt.title('Modello SIR', fontsize=20)
plt.plot(t, S, 'blue' , label='Suscettibili')
plt.plot(t, I, 'red'  , label='Infetti')
plt.plot(t, R, 'black', label='Rimessi')
plt.legend(loc='best')
plt.grid()

plt.show()