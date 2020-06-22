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

U,S,VT = m_3a.SVD()
print(U)
print(S)
print(VT)
m_3a2 = U*S
print(m_3a2)
m_3a2 *= VT
print(m_3a2)
print(m_3a)

print("Scipy version...")
y = np.asarray(m_3a)
print(y)
U2, s2, Vh = linalg.svd(y)
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

print("\n(4) Pseudo Inverse Tests...")
print("\n(4a)")

# construct an n x m matrix m_4a
m_4a = Matrix("m_1a",3,2,[1,1,1,2,3,2])
print(m_4a)

# compute the pseudo-inverse using two different approaches
# the direct way uses the function from the python bindings
print("\nDirect way")
m_4ap = m_4a.pseudo_inverse()
print(m_4ap)
# the indirect way gets the SVD decomposition from the bindings
# and then computes the matrix product.
print("\nIndirect way")
U, S, VT = m_4a.SVD()
print("\nThe matrices U, Sigma and (V)^T")
print(U)
print(S)
print(VT)
for i in range(min(S.get_num_cols(),S.get_num_rows())):
    S[i,i] = 1/S[i,i]
UT = U.T()
ST = S.T()
V = VT.T()
print("\nThe matrices (U)^T, Sigma^+ and V")
print(UT)
print(ST)
print(V)
m_4ap2 = (VT.T())*(S.T())*(U.T())
print(m_4ap2)

print("\nCheck the product")
m_4app = m_4a * (m_4ap * m_4a)
print(m_4app)

print("Scipy version...")
print("Matrix m_4a")
y = np.asarray(m_4a)
print(y)
print("Direct way")
print("Pseudo-inverse of m_4a")
yp = np.linalg.pinv(y)
print(yp)
print("Indirect way")
u, s, vh = linalg.svd(y)
a = m_4a.get_num_rows()
b = m_3a.get_num_cols()
sigma = np.zeros((a,b))
for i in range(min(a,b)):
    sigma[i, i] = s2[i]
vh = np.transpose(vh)
u = np.transpose(u)
np.transpose(sigma)
a1 = np.dot(vh, np.dot(sigma, u))
print(a1)
np.allclose(y, np.dot(y, np.dot(a1, y)))
