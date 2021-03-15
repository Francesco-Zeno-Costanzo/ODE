import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from scipy import constants as cst

step=20000
dt=1E-5
l=1E-6
I=10    #momento d'ineriza
m=10    #massa
p=0.001 #dipolo mobile
q=1E-4  #dipolo fisso
k=1/(4*np.pi*cst.value(u'vacuum electric permittivity'))


def E(x,y):
    Q_x=q*np.cos(np.pi/2)
    Q_y=q*np.sin(np.pi/2)
    Q=np.array([Q_x,Q_y])
    R=np.array([x,y])
    Ex=k*(3*(Q[0]*x+Q[1]*y)*x - (x**2+y**2)*Q[0])/(np.sqrt((x**2+y**2))**5)
    Ey=k*(3*(Q[0]*x+Q[1]*y)*y - (x**2+y**2)*Q[1])/(np.sqrt((x**2+y**2))**5)
    return Ex, Ey

def g(x, y, vx, vy, o, w):
    #............... f originaria per la forza
    P=np.array([p*np.cos(o),p*np.sin(o)])
    dEx=np.array([(E(x+l,y)[0]-E(x,y)[0])/l,(E(x+l,y)[1]-E(x,y)[1])/l])
    dEy=np.array([(E(x,y+l)[0]-E(x,y)[0])/l,(E(x,y+l)[1]-E(x,y)[1])/l])
    Fx=(P[0]*dEx[0]+P[1]*dEx[1])
    Fy=(P[0]*dEy[0]+P[1]*dEy[1])
    #...............
    x_dot=vx
    y_dot=vy
    vx_dot=Fx/m
    vy_dot=Fy/m
    #.............. f originaria pe il momento
    #P=np.array([p*np.cos(o),p*np.sin(o)])
    M=np.array([0], dtype=float)
    M[0]=(P[0]*E(x, y)[1]-P[1]*E(x, y)[0])
    #M[0]=np.cross(P,E(x,y))
    #..............
    o_dot= w
    w_dot=M[0]/I
    return x_dot, y_dot, vx_dot, vy_dot, o_dot, w_dot


xp=np.zeros(step+1)
yp=np.zeros(step+1)
vx=np.zeros(step+1)
vy=np.zeros(step+1)
wp=np.zeros(step+1)
op=np.zeros(step+1)
tp=np.zeros(step+1)
xp[0], yp[0], vx[0], vy[0], op[0], wp[0] = (1, 0, 0, 0, 0, 0)

for i in range(step):
    xk1, yk1, vxk1, vyk1, ok1, wk1 = g(xp[i], yp[i], vx[i], vy[i], op[i], wp[i])
    xk2, yk2, vxk2, vyk2, ok2, wk2 = g(xp[i]+xk1*dt/2, yp[i]+yk1*dt/2, vx[i]+vxk1*dt/2, vy[i]+vyk1*dt/2, op[i]+ok1*dt/2, wp[i]+wk1*dt/2)
    xk3, yk3, vxk3, vyk3, ok3, wk3 = g(xp[i]+xk2*dt/2, yp[i]+yk2*dt/2, vx[i]+vxk2*dt/2, vy[i]+vyk2*dt/2, op[i]+ok2*dt/2, wp[i]+wk2*dt/2)
    xk4, yk4, vxk4, vyk4, ok4, wk4 = g(xp[i]+xk3*dt, yp[i]+yk3*dt, vx[i]+vxk3*dt, vy[i]+vyk3*dt, op[i]+ok3*dt, wp[i]+wk3*dt)
    xp[i+1]=xp[i]+(dt/6)*(xk1 + 2*xk2 + 2*xk3 + xk4)
    yp[i+1]=yp[i]+(dt/6)*(yk1 + 2*yk2 + 2*yk3 + yk4)
    vx[i+1]=vx[i]+(dt/6)*(vxk1 + 2*vxk2 + 2*vxk3 + vxk4)
    vy[i+1]=vy[i]+(dt/6)*(vyk1 + 2*vyk2 + 2*vyk3 + vyk4)
    op[i+1]=op[i]+(dt/6)*(ok1 + 2*ok2 + 2*ok3 + ok4)
    wp[i+1]=wp[i]+(dt/6)*(wk1 + 2*wk2 + 2*wk3 + wk4)
    tp[i+1]=tp[i]+dt
    if xp[i+1]**2 +yp[i+1]**2 < 0.01:
        break


pE=np.zeros(step+1)
for i in range(len(xp)):
    if xp[i]**2 +yp[i]**2 > 0.01:
        pE[i]=p*np.cos(op[i])*E(xp[i],yp[i])[0]+p*np.sin(op[i])*E(xp[i],yp[i])[1]


U=m*(vx**2 + vy**2)/2 + I*(wp**2)/2 - pE

n=1000
Mx=np.max(xp)
My=np.max(yp)
if Mx and My > 15:
    xlim, ylim = 15, 15
else:
    xlim, ylim = Mx+0.5, My+0.5

a = np.linspace(-xlim,xlim,n)
b = np.linspace(-ylim,ylim,n)
x, y = np.meshgrid(a,b)

Ex, Ey = E(x,y)[0], E(x, y)[1]

u, v = Ex/np.hypot(Ex,Ey), Ey/np.hypot(Ex,Ey)

plt.figure(1)
plt.title('Interazione dipo fisso-dipolo mobile', fontsize=20)
plt.xlabel('asse x', fontsize=10)
plt.ylabel('asse y', fontsize=10)
plt.plot(xp,yp,"r", label='traiettoria')
plt.errorbar(0,0,fmt="^",markersize=6)
plt.streamplot(x, y, u, v, color="blue", density=2.5,linewidth=0.3, arrowstyle='->', arrowsize=1)
plt.legend(loc='best')


plt.figure(2)
plt.errorbar(tp, op, fmt='.', markersize=0.7, color='blue', label='angolo')
plt.errorbar(tp, wp, fmt='.', markersize=0.7, color='black', label='velocità angolare')
plt.errorbar(tp, np.sqrt(vx**2 +vy**2), fmt='.', markersize=0.7, color='red', label='modulo della velocità')
plt.legend(loc='best')
plt.grid()

plt.figure(3)
#plt.ylim(-1, 1)
plt.errorbar(tp, U-U[0], fmt='.',  label='energia totale')
plt.legend(loc='best')
plt.grid()
plt.show()
