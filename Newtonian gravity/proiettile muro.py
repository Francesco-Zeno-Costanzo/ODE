import numpy as np
import  matplotlib.pyplot  as  plt
import matplotlib.animation as animation

a=0.4
b=0.6
g=9.81
m=12.7
tf=1

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
dt = tf/num_steps

xs = np.zeros(num_steps+1)
ys = np.zeros(num_steps+1)
vy = np.zeros(num_steps+1)
vx = np.zeros(num_steps+1)
ts = np.zeros(num_steps+1)

xs[0], vx[0], ys[0], vy[0] = (0, 2, 0, 2)

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
plt.title('Proiettile contro muro di gelatina')
plt.xlabel('X(t)')
plt.ylabel('Y(t)')
plt.plot(a*p, ll, 'k')
plt.plot(b*p, ll, 'k')
plt.plot(xs, ys)
plt.grid()
plt.figure(2)
plt.plot(ts, np.sqrt(vx**2+vy**2), 'k', label='v(t)')
plt.plot(ts, vx, label='$v_x$(t)')
plt.plot(ts, vy, label='$v_y$(t)')
plt.xlabel('t')
plt.legend(loc='best')
plt.grid()



fig = plt.figure(3)

plt.xlim(np.min(xs)-0.1,np.max(xs)+0.1)
plt.ylim(np.min(ys)-0.1,np.max(ys)+0.1)

line, = plt.plot([],[],'r-')
dot, = plt.plot([], [], 'ro')

def animate(i):
    dot.set_data(xs[i],ys[i])
    line.set_data(xs[:i],ys[:i])
    return  dot, line


myAnimation = animation.FuncAnimation(fig, animate, frames=range(0, 10000, 2),interval=1, blit=True, repeat=True)

plt.grid()
plt.title('Proiettile contro muro di gelatina')
plt.xlabel('X(t)')
plt.ylabel('Y(t)')
plt.plot(a*p, ll, 'k')
plt.plot(b*p, ll, 'k')

#anim.save('pinem.mp4', fps=30, extra_args=['-vcodec', 'libx264']) #salva sulla cartella corrente l'animazione
plt.show()