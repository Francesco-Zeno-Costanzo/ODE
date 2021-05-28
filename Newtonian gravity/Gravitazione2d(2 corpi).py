import numpy as np
import matplotlib.pyplot as plt

m=100
M=m*1000
def F(x, y, x1, y1):
    r12=np.sqrt((x-x1)**2 + (y-y1)**2)
    Fx = - m*M * (x-x1)/r12**3
    Fy = - m*M * (y-y1)/r12**3
    return Fx, Fy, -Fx, -Fy

steps=40000
dt=0.0000002

X1=np.zeros((2, steps+1))
V1=np.zeros((2, steps+1))
X2=np.zeros((2, steps+1))
V2=np.zeros((2, steps+1))

X1[:, 0]=(-1, 0)
V1[:, 0]=(0, 100)
X2[:, 0]=(0, 0)
V2[:, 0]=(0, 0)

#coefficienti ruth 3
c = np.array([2/3, -2/3, 1])
d = np.array([7/24, 3/4, -1/24])

for i in range(steps):
    x, y = X1[:, i]; vx, vy = V1[:, i]
    x1, y1 = X2[:, i]; vx1, vy1 = V2[:, i]

    for k in range(len(d)):
        vx+=  (d[k]*F(x, y, x1, y1)[0]*dt)/m
        vy+=  (d[k]*F(x, y, x1, y1)[1]*dt)/m
        vx1+= (d[k]*F(x, y, x1, y1)[2]*dt)/M
        vy1+= (d[k]*F(x, y, x1, y1)[3]*dt)/M

        x+=  c[k]*vx*dt
        y+=  c[k]*vy*dt
        x1+= c[k]*vx1*dt
        y1+= c[k]*vy1*dt

    X1[:, i+1]=(x, y)
    V1[:, i+1]=(vx, vy)
    X2[:, i+1]=(x1, y1)
    V2[:, i+1]=(vx1, vy1)

x = X1[0,:]-X2[0,:]
y = X1[1,:]-X2[1,:]
U = m*(np.sum(V1**2, axis=0))/2 +M*(np.sum(V2**2, axis=0))/2 - (m*M)/np.sqrt(x**2 + y**2)

t=np.linspace(0, steps*dt, len(U))

plt.figure(1)

plt.title('Massa in campo gravitazionale', fontsize=20)
plt.xlabel('X(t)')
plt.ylabel('Y(t)')

plt.grid()
plt.plot(X1[0,:], X1[1,:], 'k')
plt.plot(X2[0,:], X2[1,:], 'r')

plt.figure(2)
plt.title('Energia del sistema: $E-E(t_0)$', fontsize=20)
plt.grid()
plt.plot(t, U-U[0])

plt.show()