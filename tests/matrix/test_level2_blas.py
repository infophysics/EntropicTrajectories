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

x = Vector('x',[1,2])
A = Matrix('A',[[1,1],[1,1]])
y = Vector('y',[1,1])
z = et.dgemv(1.1,A,x,2.0,y)
print(z)
print("\nScipy version")
x = [1,2]
A = [[1,1],[1,1]]
y = [1,1]
z = blas.dgemv(1.1,A,x,beta=2.0,y=y)
print(z)
