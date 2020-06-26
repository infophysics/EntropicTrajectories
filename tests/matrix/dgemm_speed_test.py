import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix, Vector
from scipy.linalg import blas
import time
import etraj.etraj as et

#  tests for matrix class
print("\n-------Tests for DGEMM speed vs. Scipy--------")
print("\n----------------------------------------------")

native_times = []
scipy_times = []
dims = []
N = 2
num_iters = 1000
for i in range(num_iters):

    #   generate random matrix
    A = [[np.random.normal(0,1,1)[0] for j in range(N)]
         for i in range(N)]
    B = [[np.random.normal(0,1,1)[0] for j in range(N)]
     for i in range(N)]
    A_et = Matrix('A',A)
    B_et = Matrix('B',B)

    #   run dgemm in ET
    start = time.time()
    C_et = et.dgemm(1.0,A_et,B_et)
    end = time.time()
    native_times.append(end - start)

    #   run dgemv in scipy
    start = time.time()
    C_scipy = blas.dgemm(1.0,A,B)
    end = time.time()
    scipy_times.append(end - start)

    dims.append(N)
    N += 1

#   plot the results
fig, axs = plt.subplots()
axs.plot(dims,native_times,color='k',linestyle='--',label="ET")
axs.plot(dims,scipy_times,color='r',linestyle='--',label="Scipy")
axs.set_xlabel(r"dimension ($n^2$)")
axs.set_ylabel("log time (log(sec))")
axs.set_yscale("log")
plt.legend()
plt.title(r"DGEMM: log time (log(sec)) vs. dimension ($n^2$)")
plt.savefig("DGEMM_test.png")
plt.show()
