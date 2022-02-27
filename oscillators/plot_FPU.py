import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

n = 36
m = 800000
m = m + 1

x = np.loadtxt(r'scrivi.txt')
X = np.zeros((n, m))

X = np.reshape(x, (n, m))
t = np.linspace(0, n-1, n)

'''
k = np.linspace(0, 5971, m)
for i in range(n):
    plt.plot(k, X[i,:])
plt.show()
'''
'''
T=t[len(t)-1]

sp2 = np.fft.rfft(X[:,m//2]) #fft per imput a valori reali
freq2 = np.fft.rfftfreq(len(t), T/len(t))
plt.plot(freq2, abs(sp2))
plt.show()
'''

fig = plt.figure(1)


plt.title('Problema di Fermi Pasta Ulam', fontsize=15)
plt.xlabel('distanza')
plt.ylabel('ampiezza')
plt.xlim(np.min(t)-0.1, np.max(t)+0.1)
plt.ylim(-1 - 0.1, 1 + 0.1)

dot = np.array([])
for l in range(n):
    dot = np.append(dot, plt.plot([], [], 'ro'))

def animate(i):
    for k in range(n):
        dot[k].set_data(k, X[k][i])
    return dot


anim = animation.FuncAnimation(fig, animate, frames=np.arange(0, m, 50) ,interval=1, blit=True, repeat=True)

plt.grid()

anim.save('FPU.mp4', fps=160)
plt.show()

