import numpy as np
from scipy.sparse import diags
import  matplotlib.pyplot as plt

def relax(f, y0, y1, init, N, h, max_iter):

    # second derivative matrix
    d2 = diags([1, -2, 1], [-1, 0, 1], shape=(N, N)).toarray()
    d2 = d2/h**2
    # bound values required
    yb = [y0/h**2] + [0]*int(N-2) + [y1/h**2]
    yb = np.array(yb)
    # init guess from imput
    yo = init

    for i in range(max_iter):
        df = np.zeros(d2.shape)
        s  = np.zeros(N)

        for i in range(N):
            s[i] = 1
            yr, yl = yo + h*s, yo - h*s
            df[i, :] = (f(x, yr, h, g, o02) - f(x, yl, h, g, o02) )/(2*h)
            s[:] = 0

        yn = yo - np.linalg.solve(d2-df, d2@yo - f(x, yo, h, g, o02) + yb)
        # conti distanza?
        yo = yn
    return yo



def f(t, y, h, g, o02):

    y_dot   = (y[2:] - y[:-2])/(2*h)              # second order derivative
    y_dot_0 = (- 3*y[0]  + 4*y[1]  - y[2] )/(2*h) # second order derivative on left  bound
    y_dot_n = (  3*y[-1] - 4*y[-2] + y[-3])/(2*h) # second order derivative on right bound
    # join everything together to get the second order derivative
    y_dot = np.insert(y_dot, 0,          y_dot_0)
    y_dot = np.insert(y_dot, len(y_dot), y_dot_n)

    # equation to solve
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

y = relax(f, y0, y1, y, N, h, 50)

plt.plot(x, y)

plt.show()
