import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

def lorenz(x, y, z, s=10, r=28, b=8/3):
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot

dt = 0.01
num_steps = 10000

t=np.linspace(0, num_steps*dt, num_steps+1)
def int(x0, y0, z0):
    xs = np.zeros(num_steps + 1)
    ys = np.zeros(num_steps + 1)
    zs = np.zeros(num_steps + 1)

    xs[0], ys[0], zs[0] = (x0, y0, z0)

    for i in range(num_steps):
        x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
        xs[i + 1] = xs[i] + (x_dot * dt)
        ys[i + 1] = ys[i] + (y_dot * dt)
        zs[i + 1] = zs[i] + (z_dot * dt)

    return xs, ys, zs

x, y, z = int(0., 1.0, 1.05)
x1, y1, z1 = int(0., 1., 1.06)


fig = plt.figure(1)
ax = fig.gca(projection='3d')


ax.set_xlabel("X(t)")
ax.set_ylabel("Y(t)")
ax.set_zlabel("Z(t)")
ax.set_title("Lorenz Attractor", fontsize=20)

ax.set_xlim(np.min(x),np.max(x))
ax.set_ylim(np.min(y),np.max(y))
ax.set_zlim(np.min(z),np.max(z))

line, = ax.plot([],[],'b-', lw=0.5)
dot, = ax.plot([], [], [],'ro')

line1, = ax.plot([],[],'k-', lw=0.5)
dot1, = ax.plot([], [], [],'go')
def animate(i):
    dot.set_data_3d(x[i],y[i], z[i])
    line.set_data_3d(x[:i],y[:i], z[:i])

    dot1.set_data_3d(x1[i],y1[i], z1[i])
    line1.set_data_3d(x1[:i],y1[:i], z1[:i])
    return  dot, line, dot1, line1

anim = animation.FuncAnimation(fig, animate, frames=num_steps, interval=0.1, blit=True, repeat=True)
#anim.save('Lorenz Attractor.mp4', fps=300, extra_args=['-vcodec', 'libx264'])
plt.show()
##
plt.figure(2)
R=np.sqrt(x**2 + y**2 + z**2)
R1=np.sqrt(x1**2 + y1**2 + z1**2)
d=R-R1
plt.plot(t, np.log(abs(d)))
plt.show()