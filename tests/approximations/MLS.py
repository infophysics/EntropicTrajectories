import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix, Grid, Approximator, ScalarField
import etraj.etraj as et

N = 100

x = np.random.normal(0,1,N)
y = np.random.normal(0,1,N)

data = np.vstack((x,y)).T

g = Grid("2D",2)
g.set_grid(data)

# f = cos(x)
f = np.cos(x)
f_field = ScalarField(g,f)

# get the approximator
app = f_field.get_approximator()

# generate
m = app.construct_B_matrix(g,[1,2,3],4,3)
print(m)
