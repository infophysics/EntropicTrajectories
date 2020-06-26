import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix, Vector
from scipy.linalg import blas
import etraj.etraj as et

#  tests for matrix class
print("\n-------Tests for Level 2 BLAS functions--------")
print("\n-----------------------------------------------")

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
print(">>> x = [1,2]")
print(">>> A = [[1,1],[1,1]]")
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

print("\n(2a) - DGER")
print(">>> x = Vector('x',[1,2,3])")
print(">>> y = Vector('y',[4,5])")
print(">>> alpha = 1")
print(">>> m = et.dger(alpha,x,y)")
print(">>> print(m)")
x = Vector('x',[1,2,3])
y = Vector('y',[4,5])
alpha = 1
m = et.dger(alpha,x,y)
print(m)
print("\nScipy version")
print(">>> x = [1,2,3]")
print(">>> y = [4,5]")
print(">>> alpha = 1")
print(">>> m = blas.dger(alpha,x,y)")
print(">>> print(m)")
x = [1,2,3]
y = [4,5]
alpha = 1
m = blas.dger(alpha,x,y)
print(m)

print("\n(2b) - DGER")
print(">>> x = Vector('x',[4,5])")
print(">>> y = Vector('y',[1,2,3])")
print(">>> alpha = 2")
print(">>> m = et.dger(alpha,x,y)")
print(">>> print(m)")
x = Vector('x',[4,5])
y = Vector('y',[1,2,3])
alpha = 2
m = et.dger(alpha,x,y)
print(m)
print("\nScipy version")
print(">>> x = [4,5]")
print(">>> y = [1,2,3]")
print(">>> alpha = 2")
print(">>> m = blas.dger(alpha,x,y)")
print(">>> print(m)")
x = [4,5]
y = [1,2,3]
alpha = 2
m = blas.dger(alpha,x,y)
print(m)

print("\n(2c) - DGER (in place)")
print(">>> x = Vector('x',[4,5])")
print(">>> y = Vector('y',[1,2,3])")
print(">>> m = Matrix('m',2,3,[1,2,3,4,5,6])")
print(">>> alpha = 2")
print(">>> et.dger(alpha,x,y,m)")
print(">>> print(m)")
x = Vector('x',[4,5])
y = Vector('y',[1,2,3])
m = Matrix('m',2,3,[1,2,3,4,5,6])
alpha = 2
et.dger(alpha,x,y,m)
print(m)
print("\nScipy version")
print(">>> x = [1,2,3]")
print(">>> y = [4,5]")
print(">>> m = [[1,2,3],[4,5,6]]")
print(">>> alpha = 2")
print(">>> blas.dger(alpha,x,y,a=m,overwrite_a=1)")
print(">>> print(m)")
x = [4,5]
y = [1,2,3]
m = [[1,2,3],[4,5,6]]
alpha = 2
m = blas.dger(alpha,x,y,a=m,overwrite_a=1)
print(m)
