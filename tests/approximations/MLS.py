import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix, Grid, Approximator, ScalarField
import etraj.etraj as et
from scipy import linalg
from scipy.linalg import svdvals



print("...Moving Least Squares Tests...")
print("(1) 1D MLS...")
N = 100
x = np.asarray([[np.random.normal(0,1,1)[0]] for i in range(N)])
g = Grid("1D",1)
g.set_grid(x)
# f = cos(x)
f = np.cos(x)
f_field = ScalarField(g,f)

# get the approximator
app = f_field.get_approximator()

# find nearest neighbors of the point p,
p = 4
g.query_neighbors(3)
neighbors = g.get_neighbors(p)
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
U, s2, Vh = linalg.svd(prod_np, lapack_driver='gesvd')
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


print("(2) 2D MLS...")
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
