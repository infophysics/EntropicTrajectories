import numpy as np
import matplotlib.pyplot as plt

import etraj.etraj as etraj
from etraj.etraj import ScalarField, Grid, Matrix



g = Grid("1D",1)
N = 100
xrange = [-5.0,5.0]
x = [[np.random.uniform(-5.0,5.0,1)[0]] for i in range(N)]
g.set_grid(x)
# checks for various functions

# f = cos(x)
f = np.cos(x)
f_field = ScalarField(g,f)
# try moving least squares
f_field.set_derivative("MLS")
f_p_exact = -np.sin(x)
f_p = f_field.gradient(x,f)
fig, axs = plt.subplots()
axs.plot(x,f,color='k',label='f(x) = cos(x)')
axs.plot(x,f_p,color='r',label="~f'(x)")
axs.plot(x,f_p_exact,color='m',label='df/dx = -sin(x)')
axs.set_xlabel("x")
axs.set_ylabel("f(x), f'(x)")
axs.set_title("f(x),f'(x) vs. x")
plt.show()
