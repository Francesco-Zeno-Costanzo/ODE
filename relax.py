import numpy as np
from scipy.sparse import diags
import  matplotlib.pyplot as plt



def f(t, y, h, g, o02):
    
    y_dot = (y[1:] - y[:-1])/h
    y_dot = np.insert(y_dot, len(y_dot), (y[-2]-y[-1])/h)
    y_ddot = -g*y_dot - o02*y
    return y_ddot


g   = 0.3                    # damping factor
o02 = 1                      # proper frequency squared
xi  = 0                      # left end of the interval
N   = 1000                   # number of points
y0  = 1 
xf = 5
y1 = 0.1862

x = np.linspace(xi, xf, N)
h = np.diff(x)[0]
y = (x - xi) * (y1 - y0)/(xf - xi) + y0


P = diags([1, -2, 1], [-1, 0, 1], shape=(N, N)).toarray()
P = P/h**2
yb = [y0/h**2] + [0]*int(N-2) + [y1/h**2]
plt.plot(x, y)

for i in range(10):
    df = np.zeros(P.shape)
    s  = np.zeros(N)
    
    for i in range(N):       # loop over functions
        s[i] = 1
        yr, yl = y + h*s, y - h*s
        df[i, :] = (f(x, yr, h, g, o02) - f(x, yl, h, g, o02) )/(2*h)
        s[:] = 0
        
    y -= np.linalg.solve(P-df, P@y - f(x, y, h, g, o02) +yb)
    
    plt.plot(x, y)

plt.show()
