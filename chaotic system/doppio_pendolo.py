import numpy as np
import scipy.optimize as so
import matplotlib.pyplot as plt
from matplotlib import animation

m = 1
l = 1
g = 9.81
tf = 10

def DP(theta1, theta2, ptheta1, ptheta2):

    theta1_dot = (6/(m*l**2))*(2*ptheta1 - 3*np.cos(theta1 - theta2)*ptheta2)/(16 - 9*(np.cos(theta1 - theta2))**2)
    theta2_dot = (6/(m*l**2))*(8*ptheta2 - 3*np.cos(theta1 - theta2)*ptheta1)/(16 - 9*(np.cos(theta1 - theta2))**2)

    ptheta1_dot = -((m*l**2)/2)*((theta1_dot*theta2_dot)*np.sin(theta1 - theta2) + 3*g/l*np.sin(theta1))
    ptheta2_dot = -((m*l**2)/2)*(-(theta1_dot*theta2_dot)*np.sin(theta1 - theta2) + g/l*np.sin(theta2))

    return theta1_dot, theta2_dot, ptheta1_dot, ptheta2_dot

num_steps = 10000
dt = tf/num_steps

theta1s = np.zeros(num_steps + 1)
theta2s = np.zeros(num_steps + 1)
ptheta1s = np.zeros(num_steps + 1)
ptheta2s = np.zeros(num_steps + 1)

theta1s[0], theta2s[0], ptheta1s[0], ptheta2s[0] = (np.pi/2, np.pi/2, 0, 0)

for i in range(num_steps):
    x0, x00, v0, v00 = DP(theta1s[i], theta2s[i], ptheta1s[i], ptheta2s[i])
    x1, x11, v1, v11 = DP(theta1s[i] + x0*dt/2, theta2s[i]+ x00*dt/2, ptheta1s[i] + v0*dt/2, ptheta2s[i]+ v00*dt/2)
    theta1s[i + 1] = theta1s[i] + dt*x1
    theta2s[i + 1] = theta2s[i] + dt*x11
    ptheta1s[i + 1] = ptheta1s[i] + dt*v1
    ptheta2s[i + 1] = ptheta2s[i] + dt*v11


x_1 = l*np.sin(theta1s)
y_1 = -l*np.cos(theta1s)
x_2 = x_1 + l*np.sin(theta2s)
y_2 = y_1 - l*np.cos(theta2s)


fig = plt.figure(1)
plt.title('Doppio pendolo')
plt.xlim(-2.2,2.2)
plt.ylim(-2.2,2.2)
plt.gca().set_aspect('equal', adjustable='box')

xf, yf= [0,x_1[0],x_2[0]],[0,y_1[0],y_2[0]]

line1,=plt.plot(xf, yf, linestyle='-', marker='o',color='k')
line2,=plt.plot([], [], linestyle='-',color='blue')
plt.grid()


def animate(i):
    xf[1]=x_1[i]
    yf[1]=y_1[i]
    xf[2]=x_2[i]
    yf[2]=y_2[i]
    line1.set_data(xf, yf)
    line2.set_data(x_2[:i], y_2[:i])
    return line1, line2

anim=animation.FuncAnimation(fig, animate, frames=range(0, num_steps, 3), interval=1, blit=True, repeat=True)

#anim.save('dp.mp4', fps=90, extra_args=['-vcodec', 'libx264'])

plt.show()
