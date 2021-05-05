import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def R(x, y, z, a=0.2, b=0.2, c=5.7):
    x_dot = -y -z
    y_dot = x + a*y
    z_dot = b+z*(x-c)
    return x_dot, y_dot, z_dot


dt = 0.005
num_steps = 60000

xs = np.zeros(num_steps + 1)
ys = np.zeros(num_steps + 1)
zs = np.zeros(num_steps + 1)


xs[0], ys[0], zs[0] = (0., 1., 1.05)


for i in range(num_steps):
    x_dot, y_dot, z_dot = R(xs[i], ys[i], zs[i])
    xp = xs[i] + dt*x_dot
    yp = ys[i] + dt*y_dot
    zp = zs[i] + dt*z_dot
    xc, yc, zc = R(xp, yp, zp)
    xs[i + 1] = xs[i] + (1/2)*dt*(x_dot + xc)
    ys[i + 1] = ys[i] + (1/2)*dt*(y_dot + yc)
    zs[i + 1] = zs[i] + (1/2)*dt*(z_dot + zc)


fig = plt.figure()
ax = fig.gca(projection='3d')

ax.plot(xs, ys, zs, lw=0.7)
ax.set_xlabel("X Axis")
ax.set_ylabel("Y Axis")
ax.set_zlabel("Z Axis")
ax.set_title("Rossel Attractor")

plt.show()