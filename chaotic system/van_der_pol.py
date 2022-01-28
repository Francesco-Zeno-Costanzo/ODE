import numpy as np
import  matplotlib.pyplot  as  plt

#numero di condizioni inziali per posizione e volocità
N = 10
M = 10

#tf = tempo di integrazione
tf = 14
#mu = 0.5

#range dei valori inziali
x_min = -3
x_max = -x_min
v_min = -4
v_max = -v_min

#array contenente i valori inziali
x_0 = np.linspace(x_min, x_max, N)
v_0 = np.linspace(v_min, v_max, M)

#numero di passi
num_steps = 10000


def eq(x, v, m):
    x_dot = v
    v_dot = m*(1 - x**2)*v - x
    d = np.array([x_dot, v_dot])
    return d



def RK45(x0, v0, m):

    i = 0
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

        xk1, vk1 = h[i]*eq(x[i], v[i], m)

        xk2, vk2 = h[i]*eq(x[i]+B1[1]*xk1,
                           v[i]+B1[1]*vk1, m)
        xk3, vk3 = h[i]*eq(x[i]+B1[2]*xk1 + B2[2]*xk2,
                           v[i]+B1[2]*vk1 + B2[2]*vk2, m)
        xk4, vk4 = h[i]*eq(x[i]+B1[3]*xk1 + B2[3]*xk2 + B3[3]*xk3,
                           v[i]+B1[3]*vk1 + B2[3]*vk2 + B3[3]*vk3, m)
        xk5, vk5 = h[i]*eq(x[i]+B1[4]*xk1 + B2[4]*xk2 + B3[4]*xk3 + B4[4]*xk4,
                           v[i]+B1[4]*vk1 + B2[4]*vk2 + B3[4]*vk3 + B4[4]*vk4, m)
        xk6, vk6 = h[i]*eq(x[i]+B1[5]*xk1 + B2[5]*xk2 + B3[5]*xk3 + B4[5]*xk4 + B5[5]*xk5,
                           v[i]+B1[5]*vk1 + B2[5]*vk2 + B3[5]*vk3 + B4[5]*vk4 + B5[5]*vk5, m)


        x[i+1] = x[i] + CH[0]*xk1 + CH[1]*xk2+ CH[2]*xk3+ CH[3]*xk4+ CH[4]*xk5+ CH[5]*xk6
        v[i+1] = v[i] + CH[0]*vk1 + CH[1]*vk2+ CH[2]*vk3+ CH[3]*vk4+ CH[4]*vk5+ CH[5]*vk6

        TEx = abs(CT[0]*xk1 + CT[1]*xk2+ CT[2]*xk3+ CT[3]*xk4+ CT[4]*xk5+ CT[5]*xk6)
        TEv = abs(CT[0]*vk1 + CT[1]*vk2+ CT[2]*vk3+ CT[3]*vk4+ CT[4]*vk5+ CT[5]*vk6)
        TE = np.min([TEx, TEv])

        hn = 0.9*h[i]*(tau/TE)**(1/5)
        if TE < tau:
            i += 1
        h[i] = hn

    x = x[:i]
    v = v[:i]

    return x, v


plt.figure(1, figsize=(9, 8))
plt.suptitle(f"Spazio delle fasi oscillatore di van der pol", fontsize=15)

for k in range(4):
    m = 0.5*(1 + k)
    plt.subplot(2,2,k+1)
    plt.title(f"$\mu$ = {m}")
    plt.grid()

    #risolvo per ogni condizione inizale per
    #avere l'intero spazio delle fasi
    for x0 in x_0:
        for v0 in v_0:
            a, b = RK45(x0, v0, m)
            plt.plot(a, b, 'b')

plt.show()

