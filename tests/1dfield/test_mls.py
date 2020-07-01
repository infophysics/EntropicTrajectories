import numpy as np
import matplotlib.pyplot as plt
# we'll need Matrix, UGrid, and Approximator
from etraj.etraj import Vector, Matrix, UGrid, Approximator, ScalarField, Monomial
import etraj.etraj as et
import time

# create a random one-dimensional grid
# between -pi and pi.
N = 2000
x = np.random.uniform(-np.pi,np.pi,N)

print("Creating UGrid")
g = UGrid(x)
# generate the function values for f(x) = cos(x)
f = np.cos(1.5*x)

# consider a random point in our grid
i_rand = np.random.randint(len(x))
x_rand, f_rand = x[i_rand], f[i_rand]

# select a number of nearest neighbors to use for the point
k = 6
g.query_neighbors(k)
# grab the nearest neighbors to the point x_rand
neighbors = g.get_neighbors(i_rand)
x_neighbors = [x[i] for i in neighbors]
f_neighbors = [f[i] for i in neighbors]

# let's use those neighbors to approximate the derivative
# of f(x) at x_rand.
print("Creating approximator", flush=True)
app = Approximator()
# construct the B matrix for first order in the Taylor expansion
print("Creating Taylor matrix", flush=True)
b_matrix = app.construct_taylor_matrix(g,neighbors,i_rand,3)
print(b_matrix)

# solve the least squares problem Ax = y
# where y is the f_neighbors
# and x are the coefficients
print("Creating f vector", flush=True)
v = Vector(f_neighbors)
print("Running DGELS", flush=True)
f_app = et.dgelsd(b_matrix,v)
print(f_app)

# let's approximate the derivative of the entire function
k = 10
n = 9
g.query_neighbors(k)
f_der = []
f_derr = []
f_der_true = -1.5*np.sin(1.5*x)
#y = g.get_neighbors()
print("looping...", flush=True)
for i in range(len(x)):
    print(i)
    temp_neighbors = g.get_neighbors(i)
    print(hex(id(temp_neighbors)))
    #b_matrix = app.construct_taylor_matrix(g,temp_neighbors,i,6)
    mon = Monomial(g.get_dim(),n)
    print(hex(id(g)))
    b = [mon.taylor_monomial_expansion(g[i],g[temp_neighbors[j]]) for j in range(k)]
    print(hex(id(b)))
    b_matrix = Matrix(b)
    print(hex(id(b_matrix)))
    print(b_matrix)
    print(hex(id(b_matrix)))
    x_neighbors = [x[m] for m in temp_neighbors]
    f_neighbors = [f[m] for m in temp_neighbors]
    v = Vector(f_neighbors)
    f_app = et.dgels(b_matrix,v)
    f_der.append(f_app[1])
    f_derr.append(f_app[2])
fig, axs = plt.subplots(figsize=(12,8))
axs.scatter(x,f,color='k',label='f(x) = cos(1.5*x)')
axs.scatter(x,f_der,color='r',label=r'$\approx df/dx$')
axs.scatter(x,f_derr,color='g',label=r'$\approx d^2f/dx^2$')
axs.scatter(x,f_der_true,color='m',label='-1.5*sin(1.5*x)')
axs.set_xlabel("x")
axs.set_ylabel("f(x)")
axs.set_title("f(x) vs. x")
plt.legend()
plt.show()
