import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**3-8*x**2-50*x+23

x = np.linspace(0, 5, 100)
plt.plot(x,f(x))
plt.grid()
plt.show()

def bisection (a, b, tol):
    xl= a
    xr= b
    while (np.abs(xl-xr)>= tol):
        c = (xl+xr)/2.0
        prod = f(xl)*f(c)
        if prod > tol: 
            xl = c
        else:
            if prod < tol:
                xr =c
    return c

answer = bisection (1, 4, 1e-5)
print(f'Bisection Method Gives Root at {answer}')
    