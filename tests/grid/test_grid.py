import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from etraj.etraj import Vector, Matrix, UGrid, Approximator, ScalarField, Monomial
import etraj.etraj as et
import time
import sys


print("\n-----------Tests for Vector class--------------")
print("\n-----------------------------------------------")


N = 1000
x = np.random.uniform(-np.pi,np.pi,N)
g = UGrid(x)

k = 11
n = 20
g.query_neighbors(k)
print("References to g: ", sys.getrefcount(g))
print("Location of g: ", hex(id(g)))

f = np.cos(1.5*x)
f_der = []
f_derr = []

mon = Monomial(g.get_dim(),n)
for i in range(len(x)):
    point = g[i]
    print("Point: ",point)
    neighbors = g.get_neighbors(i)
    print("Neighbors: ", neighbors)
    b = []
    for j in range(len(neighbors)):
        b.append(mon.taylor_monomial_expansion(point,g[neighbors[j]]))
    b_matrix = Matrix('b',b)
    print("B Matrix: ")
    print(b_matrix)
    f_neighbors = [f[m] for m in neighbors]
    v = Vector(f_neighbors)
    print("Vector: ")
    print(v)
    print("References to b_matrix: ", sys.getrefcount(b_matrix))
    print("Location of b_matrix: ", hex(id(b_matrix)))
    print("References to v: ", sys.getrefcount(v))
    print("Location of v: ", hex(id(v)))
    f_app = et.dgels(b_matrix,v)
    f_der.append(f_app[1])
    f_derr.append(f_app[2])

print("Location of b_matrix", hex(id(b_matrix)))
print("Location of v", hex(id(v)))
print("Location of g", hex(id(g)))
print("References to g: ", sys.getrefcount(g))
print("References to b: ", sys.getrefcount(b))
print("References to point: ", sys.getrefcount(point))
fig, axs = plt.subplots(figsize=(12,8))
axs.scatter(x,f,color='k',label='f(x) = cos(1.5*x)')
axs.scatter(x,f_der,color='r',label=r'$\approx df/dx$')
axs.scatter(x,f_derr,color='g',label=r'$\approx d^2f/dx^2$')
axs.set_xlabel("x")
axs.set_ylabel("f(x)")
axs.set_title("f(x) vs. x")
plt.legend()
plt.show()
