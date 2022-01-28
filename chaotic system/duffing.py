import numpy as np
import matplotlib.pyplot as plt

def f(t, x, y, *pars):

    delta = pars[0]
    alpha = pars[1]
    beta  = pars[2]
    gamma = pars[3]
    omega = pars[4]

    x_dot = y
    y_dot = - delta*y - alpha*x - beta*x**3 + gamma*np.sin(omega*t)

    return x_dot, y_dot


dt = 0.006
num_steps = 20000


def RK4(x0, v0, *pars):

    xs = np.zeros(num_steps + 1)
    ys = np.zeros(num_steps + 1)
    ts = np.zeros(num_steps + 1)

    xs[0], ys[0] = (x0, v0)

    for i in range(num_steps):
        xk1, yk1 = f(ts[i], xs[i], ys[i], *pars)
        xk2, yk2 = f(ts[i] + dt/2, xs[i] + xk1*dt/2, ys[i] + yk1*dt/2, *pars)
        xk3, yk3 = f(ts[i] + dt/2, xs[i] + xk2*dt/2, ys[i] + yk2*dt/2, *pars)
        xk4, yk4 = f(ts[i] + dt, xs[i] + xk3*dt, ys[i] + yk3*dt, *pars)
        xs[i + 1] = xs[i] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
        ys[i + 1] = ys[i] + (dt/6)*(yk1 + 2*yk2 + 2*yk3 + yk4)
        ts[i + 1] = ts[i] + dt

    return xs, ys, ts

alpha, beta, delta, omega = -1, 1, 0.3, 1.2
plt.figure(1, figsize=(8, 10))
plt.suptitle('Subharmonics of the Duffing equation with: \n '
fr'$\alpha$ = {alpha}, $\beta$ = {beta}, $\delta$ = {delta}, $\omega$ = {omega}')
G = np.array([0.20, 0.29, 0.37, 0.50, 0.65])
i = 1
for g in G:
    x, y, t = RK4(1, 0, delta, alpha, beta, g, omega)
    plt.subplot(len(G), 2, i)
    plt.ylabel(f'$\gamma$ = {g}')
    plt.grid()
    plt.plot(t, x)
    plt.subplot(len(G), 2, i+1)
    plt.grid()
    plt.plot(x, y)
    i += 2

plt.show()