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
