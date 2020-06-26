import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Vector
from scipy.linalg import blas
import time
import etraj.etraj as et

#  tests for matrix class
print("\n-------Tests for DDOT speed vs. Scipy--------")
print("\n----------------------------------------------")

native_times = []
scipy_times = []
dims = []
N = 1
num_iters = 1000
for i in range(num_iters):

    #   generate random vector
    x_vec = np.random.normal(0,1,N)
    x_et = Vector('x',x_vec)

    y_vec = np.random.normal(0,1,N)
    y_et = Vector('y',y_vec)

    #   run dgemm in ET
    start = time.time()
    s = et.ddot(x_et,y_et)
    end = time.time()
    native_times.append(end - start)

    #   run dgemv in scipy
    start = time.time()
    s = blas.ddot(x_vec,y_vec)
    end = time.time()
    scipy_times.append(end - start)

    dims.append(N)
    N += 100

#   plot the results
fig, axs = plt.subplots()
axs.plot(dims,native_times,color='k',linestyle='--',label="ET")
axs.plot(dims,scipy_times,color='r',linestyle='--',label="Scipy")
axs.set_xlabel(r"dimension ($n$)")
axs.set_ylabel("log time (log(sec))")
axs.set_yscale("log")
plt.legend()
plt.title(r"DDOT: log time (log(sec)) vs. dimension ($n$)")
plt.savefig("DDOT_test.png")
plt.show()
