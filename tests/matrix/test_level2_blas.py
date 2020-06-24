import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix, Vector
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
