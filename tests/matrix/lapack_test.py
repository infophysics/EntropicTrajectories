import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu, lu_factor
from etraj.etraj import Matrix


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
