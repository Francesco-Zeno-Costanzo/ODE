import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D

def L(x, y, z, s=10, r=28, b=8/3):
    x_dot = s*(y - x)
    y_dot = r*x - y - x*z
    z_dot = x*y - b*z
    return x_dot, y_dot, z_dot

dt = 0.005
num_steps = 10000

t=np.linspace(0, num_steps*dt, num_steps+1)

def int(x0, y0, z0):
    xs = np.zeros(num_steps + 1)
    ys = np.zeros(num_steps + 1)
    zs = np.zeros(num_steps + 1)

    xs[0], ys[0], zs[0] = (x0, y0, z0)

    for i in range(num_steps):
        x_dot, y_dot, z_dot = L(xs[i], ys[i], zs[i])
        xp = xs[i] + dt*x_dot
        yp = ys[i] + dt*y_dot
        zp = zs[i] + dt*z_dot
        xc, yc, zc = L(xp, yp, zp)
        xs[i + 1] = xs[i] + (1/2)*dt*(x_dot + xc)
        ys[i + 1] = ys[i] + (1/2)*dt*(y_dot + yc)
        zs[i + 1] = zs[i] + (1/2)*dt*(z_dot + zc)

    return xs, ys, zs

x, y, z = int(0., 1.0, 1.05)
x1, y1, z1 = int(0., 1., 1.06)


plt.figure(1)
plt.suptitle('singole traiettorie')
plt.subplot(311)
plt.ylabel('x(t)')
plt.xlabel('t')
plt.plot(t, x)
plt.grid()
plt.subplot(312)
plt.ylabel('y(t)')
plt.xlabel('t')
plt.plot(t, y)
plt.grid()
plt.subplot(313)
plt.ylabel('z(t)')
plt.xlabel('t')
plt.plot(t, z)
plt.grid()

fig = plt.figure(2)
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

anim = animation.FuncAnimation(fig, animate, frames=range(0, num_steps, 2), interval=1, blit=True, repeat=True)
#anim.save('Lorenz Attractor.mp4', fps=300, extra_args=['-vcodec', 'libx264'])
plt.show()
##

R0 = np.sqrt(x**2 + y**2 + z**2)
R1 = np.sqrt(x1**2 + y1**2 + z1**2)
dist = R0 - R1

plt.figure(3)
plt.title('Differenza fra due soluzioni')
plt.yscale('log')
plt.plot(t, abs(dist))
plt.grid()
plt.show()