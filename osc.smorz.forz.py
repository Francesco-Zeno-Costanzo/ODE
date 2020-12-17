import numpy as np
import matplotlib.pyplot as plt

def f(t, x, y, b=0.1, c=1, d=1, A=1):
    x_dot = y
    y_dot = -b*y - c*x + A*np.sin(d*t)
    return x_dot, y_dot


dt = 0.005
num_steps = 30000


xs = np.zeros(num_steps + 1)
ys = np.zeros(num_steps + 1)
ts = np.zeros(num_steps + 1)

xs[0], ys[0]=(0, 0)


for i in range(num_steps):
    xk1, yk1 = f(ts[i], xs[i], ys[i])
    xk2, yk2 = f(ts[i]+dt/2, xs[i] + xk1*dt/2, ys[i]+ yk1*dt/2)
    xk3, yk3 = f(ts[i]+dt/2, xs[i] + xk2*dt/2, ys[i]+ yk2*dt/2)
    xk4, yk4 = f(ts[i]+dt, xs[i] + xk3*dt, ys[i]+ yk3*dt)
    xs[i + 1] = xs[i] + (dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
    ys[i + 1] = ys[i] + (dt/6)*(yk1 + 2*yk2 + 2*yk3 + yk4)
    ts[i + 1] = ts[i] + dt


plt.figure(1)
plt.title('Oscillatore smorzato forzato', fontsize=20)
plt.ylabel('Ampiezza', fontsize=10)
plt.xlabel('tempo', fontsize=10)

plt.plot(ts, xs, 'blue')

plt.grid()
plt.show()