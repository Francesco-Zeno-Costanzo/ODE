import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

m=100
M=m*1000
def F(x, y, z, x1, y1, z1):
    r12=np.sqrt((x-x1)**2 + (y-y1)**2 + (z-z1)**2)
    Fx = - m*M * (x-x1)/r12**3
    Fy = - m*M * (y-y1)/r12**3
    Fz = - m*M * (z-z1)/r12**3
    return Fx, Fy, Fz, -Fx, -Fy, -Fz

steps=20000
dt=0.000001

X1=np.zeros((3, steps+1))
V1=np.zeros((3, steps+1))
X2=np.zeros((3, steps+1))
V2=np.zeros((3, steps+1))

X1[:, 0]=(-1, 0, 1)
V1[:, 0]=(0, 100, 0)
X2[:, 0]=(0, 0, 0)
V2[:, 0]=(0, 0, 0)

#coefficienti ruth 3
c = np.array([2/3, -2/3, 1])
d = np.array([7/24, 3/4, -1/24])

for i in range(steps):
    x, y, z = X1[:, i]; vx, vy, vz = V1[:, i]
    x1, y1, z1 = X2[:, i]; vx1, vy1, vz1 = V2[:, i]

    for k in range(len(d)):
        vx+=  (d[k]*F(x, y, z, x1, y1, z1)[0]*dt)/m
        vy+=  (d[k]*F(x, y, z, x1, y1, z1)[1]*dt)/m
        vz+=  (d[k]*F(x, y, z, x1, y1, z1)[2]*dt)/m
        vx1+= (d[k]*F(x, y, z, x1, y1, z1)[3]*dt)/M
        vy1+= (d[k]*F(x, y, z, x1, y1, z1)[4]*dt)/M
        vz1+= (d[k]*F(x, y, z, x1, y1, z1)[5]*dt)/M

        x+=  c[k]*vx*dt
        y+=  c[k]*vy*dt
        z+=  c[k]*vz*dt
        x1+= c[k]*vx1*dt
        y1+= c[k]*vy1*dt
        z1+= c[k]*vz1*dt

    X1[:, i+1]=(x, y, z)
    V1[:, i+1]=(vx, vy, vz)
    X2[:, i+1]=(x1, y1, z1)
    V2[:, i+1]=(vx1, vy1, vz1)

x = X1[0,:]-X2[0,:]
y = X1[1,:]-X2[1,:]
z = X1[2,:]-X2[2,:]
U = m*(np.sum(V1**2, axis=0))/2 +M*(np.sum(V2**2, axis=0))/2 - (m*M)/np.sqrt(x**2 + y**2 + z**2)

t=np.linspace(0, steps*dt, len(U))

plt.figure(1)
plt.title('Energia del sistema: $E-E(t_0)$', fontsize=20)
plt.grid()
plt.plot(t, U-U[0])


fig1 = plt.figure(2)
ax = fig1.gca(projection='3d')
plt.grid()

ax.set_title('Massa in campo gravitazionale', fontsize=20)
ax.set_xlabel('X(t)')
ax.set_ylabel('Y(t)')
ax.set_zlabel('Z(t)')

ax.plot(X1[0,:], X1[1,:], X1[2,:], lw=3)
ax.plot(X2[0,:], X2[1,:], X2[2,:], lw=3)

plt.show()