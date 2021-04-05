import numpy as np
import  matplotlib.pyplot  as  plt

a=0.4
b=0.6
g=9.81
m=4.7
xf=1.5

def F(vx, vy, x):
    if x < a:
        Fy = -g
        Fx = 0
    elif x>a and x<b:
        Fy = -m*vy -g
        Fx = -m*vx
    elif x>b:
        Fy = -g
        Fx = 0
    return Fy, Fx

c = np.array([2/3, -2/3, 1])
d = np.array([7/24, 3/4, -1/24])

num_steps = 10000
dt = xf/num_steps

xs = np.zeros(num_steps+1)
ys = np.zeros(num_steps+1)
vy = np.zeros(num_steps+1)
vx = np.zeros(num_steps+1)
ts = np.zeros(num_steps+1)

xs[0], vx[0], ys[0], vy[0] = (0, 1, 10, 0)

for i in range(num_steps):
    x, v1x, y, v1y = xs[i], vx[i], ys[i], vy[i]
    for k in range(len(d)):
        v1x+=d[k]*F(v1x, v1y, x)[1]*dt
        v1y+=d[k]*F(v1x, v1y, x)[0]*dt
        x+=c[k]*v1x*dt
        y+=c[k]*v1y*dt
    xs[i + 1] = x
    ys[i + 1] = y
    vx[i + 1] = v1x
    vy[i + 1] = v1y
    ts[i + 1] = ts[i] + dt

p=np.ones(1000)
ll=np.linspace(np.min(ys)-0.1, np.max(ys)+0.1, 1000)

plt.figure(1)
plt.title('proiettile contro muro di gelatina')
plt.plot(a*p, ll, 'k')
plt.plot(b*p, ll, 'k')
plt.plot(xs, ys)
plt.grid()
plt.figure(2)
plt.plot(ts, np.sqrt(vx**2+vy**2), 'k', label='v')
plt.plot(ts, vx, label='$v_x$')
plt.plot(ts, vy, label='$v_y$')
plt.legend(loc='best')
plt.grid()
plt.show()