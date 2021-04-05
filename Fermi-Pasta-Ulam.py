import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

a=1.

def F(x):
    F=np.zeros(len(x))
    F[0]=0
    for k in range(1, len(x)-1):
        F[k]=(-2*x[k]+x[k-1]+x[k+1])*(1+a*(x[k+1]-x[k-1]))
    F[len(x)-1]= 0
    return F

xf=391.455
num_steps = 100000
dt = xf/num_steps
N=13
x = np.zeros((num_steps+1, N))
v = np.zeros((num_steps+1, N))

ts3 = np.zeros(num_steps+1)
for j in range(N):
    x[0, j]=np.sin(0.261*j)
    v[0, j]=0


for i in range(num_steps):
    F0 = F(x[i][:])
    for j in range(N):
        x[i + 1][j] = x[i][j] + v[i][j]*dt + F0[j]*((dt**2)/2)
    F1 = F(x[i+1][:])
    for j in range(N):
        v[i + 1][j] = v[i][j] + (F0[j]+F1[j])*dt/2
    ts3[i + 1] = ts3[i] + dt

t=np.linspace(0, N-1, N)
fig = plt.figure()
'''
for l in range(N):
    plt.plot(ts3, x[:,l])
'''
plt.title('Problema di Fermi Pasta Ulam', fontsize=15)
plt.xlabel('distanza')
plt.ylabel('ampiezza')
plt.xlim(np.min(t)-0.1, np.max(t)+0.1)
plt.ylim(np.min(x)-0.1, np.max(x)+0.1)
dot=np.array([])
for l in range(N):
    dot=np.append(dot, plt.plot([], [], 'ro'))

def animate(i):
    for k in range(N):
        dot[k].set_data(k, x[i][k])
    return dot[0], dot[1], dot[2], dot[3], dot[4], dot[5], dot[6], dot[7], dot[8], dot[9], dot[10], dot[11], dot[12]

#anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, num_steps+1, 10) ,interval=1, blit=True, repeat=True)

plt.grid()

anim.save('FPU.mp4', fps=160,  extra_args=['-vcodec', 'libx264'])
plt.show()

