import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix, Vector
from scipy.linalg import blas
import time
import etraj.etraj as et

#  tests for matrix class
print("\n-------Tests for DGEMV speed vs. Scipy--------")
print("\n----------------------------------------------")

native_times = []
scipy_times = []
dims = []
N = 2
num_iters = 110
for i in range(num_iters):
    #   generate random vector
    x_vec = np.random.normal(0,1,N)
    x_et = Vector('x',x_vec)
    #   generate random matrix
    A = [[np.random.normal(0,1,1)[0] for j in range(N)]
         for i in range(N)]
    A_et = Matrix('A',A)
    alpha = 1.0

    #   run dgemv in ET
    start = time.time()
    y_et = et.dgemv(alpha,A_et,x_et)
    end = time.time()
    native_times.append(end - start)

    #   run dgemv in scipy
    start = time.time()
    y_scipy = blas.dgemv(alpha,A,x_vec)
    end = time.time()
    scipy_times.append(end - start)

    dims.append(N)
    N += 1

#   plot the results
fig, axs = plt.subplots()
axs.plot(dims,native_times,color='k',linestyle='--',label="ET")
axs.plot(dims,scipy_times,color='r',linestyle='--',label="Scipy")
axs.set_xlabel("Dimension")
axs.set_ylabel("Time (sec)")
plt.legend()
plt.title("DGEMV: Time (sec) vs. Dimension")
plt.show()
