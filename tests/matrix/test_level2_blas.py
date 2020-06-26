import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix, Vector
from scipy.linalg import blas
import etraj.etraj as et

#  tests for matrix class
print("\n-------Tests for Level 2 BLAS functions--------")
print("\n-----------------------------------------------")

print("\n(1) xGEMV")
print("\n(1a) - DGEMV")
print(">>> x = Vector('x',[1,2])")
print(">>> A = Matrix('A',[[1,1],[1,1]])")
print(">>> alpha = 1.1")
print(">>> y = et.dgemv(alpha,A,x)")
print(">>> print(y)")
x = Vector('x',[1,2])
A = Matrix('A',[[1,1],[1,1]])
alpha = 1.1
y = et.dgemv(alpha,A,x)
print(y)
print("\nScipy version")
print(">>> y_scipy = blas.dgemv(alpha,A,x)")
print(">>> print(y_scipy)")
x = [1,2]
A = [[1,1],[1,1]]
y_scipy = blas.dgemv(alpha,A,x)
print(y_scipy)

print("\n(1b) - DGEMV with y and beta")
print(">>> x = Vector('x',[1,2])")
print(">>> A = Matrix('A',[[1,1],[1,1]])")
print(">>> alpha = 1.1")
print(">>> beta = 2.0")
print(">>> y = Vector('y',[2,1])")
print(">>> et.dgemv(alpha,A,x,beta,y)")
print(">>> print(y)")
x = Vector('x',[1,2])
A = Matrix('A',[[1,1],[1,1]])
alpha = 1.1
beta = 2.0
y = Vector('y',[2,1])
et.dgemv(alpha,A,x,beta,y)
print(y)
print("\nScipy version")
print(">>> x = [1,2]")
print(">>> A = [[1,1],[1,1]]")
print(">>> y = [2,1]")
print(">>> y = blas.dgemv(alpha,A,x,beta=beta,y=y,overwrite_y=1)")
print(">>> print(y)")
x = [1,2]
A = [[1,1],[1,1]]
y = [2,1]
y = blas.dgemv(alpha,A,x,beta=beta,y=y,overwrite_y=1)
print(y)
