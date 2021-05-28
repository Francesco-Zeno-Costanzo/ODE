import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
m1=100
m=100
M=100000

def F(x, y, z, x1, y1, z1, x2, y2, z2):
    r12=np.sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)
    r13=np.sqrt((x-x2)**2 + (y-y2)**2 + (z-z2)**2)
    r23=np.sqrt((x1-x2)**2 + (y1-y2)**2+ (z1-z2)**2)

    Fx= - m*M * (x-x1)/r12**3 - m*m1 * (x-x2)/r13**3
    Fy= - m*M * (y-y1)/r12**3 - m*m1 * (y-y2)/r13**3
    Fz= - m*M * (z-z1)/r12**3 - m*m1 * (z-z2)/r13**3

    F1x= - m*M * (x1-x)/r12**3 - m1*M * (x1-x2)/r23**3
    F1y= - m*M * (y1-y)/r12**3 - m1*M * (y1-y2)/r23**3
    F1z= - m*M * (z1-z)/r12**3 - m1*M * (z1-z2)/r23**3

    F2x= - m1*M * (x2-x1)/r23**3 - m*m1 * (x2-x)/r13**3
    F2y= - m1*M * (y2-y1)/r23**3 - m*m1 * (y2-y)/r13**3
    F2z= - m1*M * (z2-z1)/r23**3 - m*m1 * (z2-z)/r13**3

    return Fx, Fy, Fz, F1x, F1y, F1z, F2x, F2y, F2z

N=3

steps=12000
dt=0.000005

X=np.zeros((3, steps+1, N))
V=np.zeros((3, steps+1, N))

#massa m
X[:, 0, 0]=(-1, 0, 1)
V[:, 0, 0]=(0, 240, 0)
#massa M
X[:, 0, 1]=(0, 0, 0)
V[:, 0, 1]=(0 ,0, 0)
#massa m1
X[:, 0, 2]=(-1, 0, -1)
V[:, 0, 2]=(0, 200, 0)
#coefficienti ruth 3
c = np.array([2/3, -2/3, 1])
d = np.array([7/24, 3/4, -1/24])

for i in range(steps):
    x,  y,  z  = X[:, i, 0]; vx,  vy,  vz  = V[:, i, 0]
    x1, y1, z1 = X[:, i, 1]; vx1, vy1, vz1 = V[:, i, 1]
    x2, y2, z2 = X[:, i, 2]; vx2, vy2, vz2 = V[:, i, 2]

    for k in range(len(d)):
        vx+=  (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[0]*dt)/m
        vy+=  (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[1]*dt)/m
        vz+=  (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[2]*dt)/m
        vx1+= (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[3]*dt)/M
        vy1+= (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[4]*dt)/M
        vz1+= (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[5]*dt)/M
        vx2+= (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[6]*dt)/m1
        vy2+= (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[7]*dt)/m1
        vz2+= (d[k]*F(x, y, z, x1, y1, z1, x2, y2, z2)[8]*dt)/m1

        x+=  c[k]*vx*dt; x1+= c[k]*vx1*dt; x2+= c[k]*vx2*dt
        y+=  c[k]*vy*dt; y1+= c[k]*vy1*dt; y2+= c[k]*vy2*dt
        z+=  c[k]*vz*dt; z1+= c[k]*vz1*dt; z2+= c[k]*vz2*dt

    X[:, i+1, 0]=(x, y, z)
    V[:, i+1, 0]=(vx, vy, vz)
    X[:, i+1, 1]=(x1, y1, z1)
    V[:, i+1, 1]=(vx1, vy1, vz1)
    X[:, i+1, 2]=(x2, y2, z2)
    V[:, i+1, 2]=(vx2, vy2, vz2)

x = X[0, :, 0] - X[0, :, 1]; x1 = X[0, :, 0] - X[0, :, 2]; x2 = X[0, :, 1] - X[0, :, 2]
y = X[1, :, 0] - X[1, :, 1]; y1 = X[1, :, 0] - X[1, :, 2]; y2 = X[1, :, 1] - X[1, :, 2]
z = X[2, :, 0] - X[2, :, 1]; z1 = X[2, :, 0] - X[2, :, 2]; z2 = X[2, :, 1] - X[2, :, 2]


T = m*(np.sum(V[:, :, 0]**2, axis=0))/2 + M*(np.sum(V[:, :, 1]**2, axis=0))/2 + m1*(np.sum(V[:, :, 2]**2, axis=0))/2
U = T - (m*M)/np.sqrt(x**2+y**2+ z**2) - (m1*m)/np.sqrt(x1**2+y1**2+z1**2) - (m1*M)/np.sqrt(x2**2+y2**2+z1**2)

t=np.linspace(0, steps*dt, len(U))

plt.figure(1)
plt.title('Energia del sistema: $E-E(t_0)$', fontsize=20)
plt.grid()
plt.plot(t, U-U[0])


fig = plt.figure(2)
ax = fig.gca(projection='3d')
plt.grid()

plt.title('Sistema a tre corpi', fontsize=20)
ax.set_xlabel('X(t)')
ax.set_ylabel('Y(t)')
ax.set_zlabel('Z(t)')

plt.grid()
colors = plt.cm.jet(np.linspace(0, 1, 3))
for l, c in zip(range(N), colors):
    ax.plot(X[0,:,l], X[1,:,l], X[2,:,l], c=c)


fig1 = plt.figure(3)
ax1 = fig1.gca(projection='3d')
plt.grid()
ax1.set_xlim(np.min(X[0,:,:]), np.max(X[0,:,:]))
ax1.set_ylim(np.min(X[1,:,:]), np.max(X[1,:,:]))
ax1.set_zlim(np.min(X[2,:,:]), np.max(X[2,:,:]))
colors = plt.cm.jet(np.linspace(0, 1, N))

dot=np.array([])
line=np.array([])

for l, c in zip(range(N), colors):
    dot=np.append(dot, plt.plot([], [], [], 'o', c=c))
    line=np.append(line, plt.plot([], [], [], '-', c=c))

def animate(i):
    for k in range(3):
        if i>5000:
            line[k].set_data_3d(X[0, i-5000:i, k], X[1, i-5000:i, k], X[2, i-5000:i, k])
        else:
            line[k].set_data_3d(X[0, :i, k], X[1, :i, k], X[2, :i, k])

        dot[k].set_data_3d(X[0, i, k], X[1, i, k], X[2, i, k])
    return dot[0], dot[1], dot[2], line[0], line[1], line[2]



anim = animation.FuncAnimation(fig, animate, frames=(steps+1), interval=1, blit=True, repeat=True)

plt.title('Sistema a tre corpi', fontsize=20)
ax1.set_xlabel('X(t)')
ax1.set_ylabel('Y(t)')
ax1.set_zlabel('Z(t)')

#anim.save('grav1.mp4', fps=120, extra_args=['-vcodec', 'libx264']) #salva sulla cartella corrente l'animazione

plt.show()