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


xs = np.zeros(num_steps + 1)
ys = np.zeros(num_steps + 1)
zs = np.zeros(num_steps + 1)

# valori inizali
xs[0], ys[0], zs[0] = (0., 1., 1.05)

#discretizzazione
for i in range(num_steps):
    x_dot, y_dot, z_dot = lorenz(xs[i], ys[i], zs[i])
    xs[i + 1] = xs[i] + (x_dot * dt)
    ys[i + 1] = ys[i] + (y_dot * dt)
    zs[i + 1] = zs[i] + (z_dot * dt)


fig = plt.figure()
ax = fig.gca(projection='3d')


ax.set_xlabel("X(t)")
ax.set_ylabel("Y(t)")
ax.set_zlabel("Z(t)")
ax.set_title("Lorenz Attractor", fontsize=20)

ax.set_xlim(np.min(xs),np.max(xs))
ax.set_ylim(np.min(ys),np.max(ys))
ax.set_zlim(np.min(zs),np.max(zs))

line, = plt.plot([],[],'b-', lw=0.5)
dot, = ax.plot([], [], [],'ro')
def animate(i):
    dot.set_data_3d(xs[i],ys[i], zs[i])
    line.set_data_3d(xs[:i],ys[:i], zs[:i])
    return  dot, line

anim = animation.FuncAnimation(fig, animate, frames=num_steps, interval=0.1, blit=True, repeat=True)
#anim.save('Lorenz Attractor.mp4', fps=300, extra_args=['-vcodec', 'libx264'])

fig1 = plt.figure()
ax1 = fig1.gca(projection='3d')
plt.grid()

ax1.set_title("Lorenz Attractor", fontsize=20)
ax1.set_xlabel('X(t)')
ax1.set_ylabel('Y(t)')
ax1.set_zlabel('Z(t)')
plt.grid()

ax1.plot(xs, ys, zs, lw=0.7)

plt.show()