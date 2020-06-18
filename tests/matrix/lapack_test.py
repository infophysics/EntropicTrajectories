import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import lu
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
m_2a = Matrix("m_2a",2,2,[-1,-2,1,2])
perm = m_2a.LU();
x = np.asarray(m_2a)
print(m_2a)
print(x)
p, l, u = lu(x)
print(perm)
print(p)
