import numpy as np
import matplotlib.pyplot as plt


o0 = 9
v0 = 0
x0 = 1

def g(x, v):
    x_dot = v
    v_dot = - o0*x
    d = np.array([x_dot, v_dot])
    return d

tf = 10
num_steps = 10000
i=0
x = np.zeros(num_steps+1)
v = np.zeros(num_steps+1)
t = np.zeros(num_steps+1)
h = np.zeros(num_steps+1)

x[0], v[0], h[0] = x0, v0, tf/num_steps

tau = 1e-10

B1 = np.array([0, 2/9, 1/12, 69/128, -17/12, 65/432])
B2 = np.array([0, 0, 1/4, -243/128, 27/4, -5/16])
B3 = np.array([0, 0, 0, 135/64, -27/5, 13/16])
B4 = np.array([0, 0, 0, 0, 16/15, 4/27])
B5 = np.array([0, 0, 0, 0, 0, 5/144])
CH = np.array([47/450, 0, 12/25, 32/225, 1/30, 6/25])
CT = np.array([-1/150, 0, 3/100, -16/75, -1/20, 6/25])

while t[i]<tf:

    t[i+1] = t[i] + h[i]

    xk1, vk1 = h[i]*g(x[i], v[i])

    xk2, vk2 = h[i]*g(x[i]+B1[1]*xk1,
                      v[i]+B1[1]*vk1)
    xk3, vk3 = h[i]*g(x[i]+B1[2]*xk1 + B2[2]*xk2,
                      v[i]+B1[2]*vk1 + B2[2]*vk2)
    xk4, vk4 = h[i]*g(x[i]+B1[3]*xk1 + B2[3]*xk2 + B3[3]*xk3,
                      v[i]+B1[3]*vk1 + B2[3]*vk2 + B3[3]*vk3)
    xk5, vk5 = h[i]*g(x[i]+B1[4]*xk1 + B2[4]*xk2 + B3[4]*xk3 + B4[4]*xk4,
                      v[i]+B1[4]*vk1 + B2[4]*vk2 + B3[4]*vk3 + B4[4]*vk4)
    xk6, vk6 = h[i]*g(x[i]+B1[5]*xk1 + B2[5]*xk2 + B3[5]*xk3 + B4[5]*xk4 + B5[5]*xk5,
                      v[i]+B1[5]*vk1 + B2[5]*vk2 + B3[5]*vk3 + B4[5]*vk4 + B5[5]*vk5)

    x[i+1] = x[i] + CH[0]*xk1 + CH[1]*xk2+ CH[2]*xk3+ CH[3]*xk4+ CH[4]*xk5+ CH[5]*xk6
    v[i+1] = v[i] + CH[0]*vk1 + CH[1]*vk2+ CH[2]*vk3+ CH[3]*vk4+ CH[4]*vk5+ CH[5]*vk6

    TE = abs(CT[0]*xk1 + CT[1]*xk2+ CT[2]*xk3+ CT[3]*xk4+ CT[4]*xk5+ CT[5]*xk6)

    hn = 0.9*h[i]*(tau/TE)**(1/5)

    if TE < tau:
        i += 1

    h[i] = hn

def Sol(t):
    return v0/np.sqrt(o0) * np.sin(np.sqrt(o0)*t) + x0*np.cos(np.sqrt(o0)*t)

def U(v, x):
    return (v**2 + o0*x**2)-(v[0]**2 + o0*x[0]**2)

x = x[:i]
v = v[:i]
t = t[:i]
h = h[:i]

plt.figure(1)
plt.subplot(221)
plt.title("Soluzione con il metodo RK4(5)", fontsize=15)
plt.grid()
plt.plot(t, x)
plt.subplot(222)
plt.title("Andamento del passo di integrazione", fontsize=15)
plt.plot(t, h)
plt.grid()
plt.subplot(223)
plt.title("Differenza tra soluzione esatta e numerica", fontsize=15)
plt.grid()
plt.plot(t, x-Sol(t))
plt.subplot(224)
plt.title("Differenza fra enerigia iniziale  ed energia al tempo t ", fontsize=15)
plt.grid()
plt.plot(t, U(v, x))
plt.show()