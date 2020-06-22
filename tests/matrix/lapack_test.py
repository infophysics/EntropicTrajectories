import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu, lu_factor
from etraj.etraj import Matrix
from scipy import linalg
from scipy.linalg import svdvals

print("...LAPACK tests...")

print("\n(1) Inverse tests...")
print("\n(1a) - m_1a_i = m_1a.inverse()")
m_1a = Matrix("m_1a",2,3,[1,1,1,2,3,2])
m_1a_i = m_1a.inverse()
print(m_1a)
print(m_1a_i)
print(m_1a_i*m_1a)


print("\n(2) LU tests...")
print("\n(2a) - permutation matrix")
m_2a = Matrix("m_2a",3,3,[1,2,3,4,5,6,7,8,9])
print(m_2a)
perm, L, U = m_2a.LU();
x = np.asarray(m_2a)
print("Numpy version")
print(x)
print("Performing PLU decomposition")
p, l, u = lu(x)
print("\nPermutation matrices")
print(perm)
print("Numpy version")
print(p)
print("\nL matrix")
print(L)
print("Numpy version")
print(l)
print("\nU Matrix")
print(U)
print("Numpy version")
print(u)


print("\n(3) SVD tests...")
n = 3
m = 3
x = [[np.random.normal(0,1,1)[0] for j in range(m)] for i in range(n)]
m_3a = Matrix("m_3a",x)
print(m_3a)
y = np.asarray(m_3a)
print(y)
U,S,VT = m_3a.SVD()
U2, s2, Vh = linalg.svd(y)
print(U)
print(S)
print(VT)
m_3a2 = U*S
print(m_3a2)
m_3a2 *= VT
print(m_3a2)
print(m_3a)
print("Scipy version...")
print("U matrix")
print(U2)
a = m_3a.get_num_rows()
b = m_3a.get_num_cols()
sigma = np.zeros((a,b))
for i in range(min(a,b)):
    sigma[i, i] = s2[i]
print("Sigma matrix")
print(sigma)
print("VT matrix")
print(Vh)
print(y)
a1 = np.dot(U2, np.dot(sigma, Vh))
print(a1)
