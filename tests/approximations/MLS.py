import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix, Grid, Approximator, ScalarField
import etraj.etraj as et
from scipy import linalg
from scipy.linalg import svdvals

N = 100

x = np.random.normal(0,1,N)
y = np.random.normal(0,1,N)
#x = np.arange(0,N,1)
#y = np.arange(0,N,1)

data = np.vstack((x,y)).T

g = Grid("2D",2)
g.set_grid(data)

# f = cos(x)
f = np.cos(x)
f_field = ScalarField(g,f)

# get the approximator
app = f_field.get_approximator()

# find nearest neighbors of the point p,
p = 4
g.query_neighbors(3)
neighbors = g.get_neighbors(p)
print(neighbors)

# generate
m = app.construct_B_matrix(g,neighbors,p,2)
print(m)
m_T = m.T()
print(m_T)
prod = m_T * m
print(prod)
s = prod.get_singular_values()

#   check singular values
a = prod.get_num_rows()
b = prod.get_num_cols()
prod_np = np.asarray(prod)
U1, S1, VT1 = prod.SVD()
U, s2, Vh = linalg.svd(prod_np)
print(U1)
print(S1)
print(VT1)
x = U1 * S1 * VT1
print(x)
print("Scipy version...")
print("U matrix")
print(U)
sigma = np.zeros((a,b))
for i in range(min(a,b)):
    sigma[i, i] = s2[i]
print("Sigma matrix")
print(sigma)
print("VT matrix")
print(Vh)
a1 = np.dot(U, np.dot(sigma, Vh))
print(x)
print(a1)
r = np.linalg.matrix_rank(a1)
print(r)
# prod_i = prod.inverse()
# print(prod_i)
# prod_i_m_T = prod_i * m_T
# print(prod_i_m_T)
#
# # play with k and n
# m2 = app.construct_B_matrix(g,[1,2,3,5,6],4,2)
# print(m2)
