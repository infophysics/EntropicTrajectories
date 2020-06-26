import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix, Vector
from scipy.linalg import blas
import etraj.etraj as et

#  tests for matrix class
print("\n-------Tests for Level 3 BLAS functions--------")
print("\n-----------------------------------------------")

print("\n(1a) - DGEMM")
print(">>> m1 = Matrix('m1',[[1,2],[3,4]])")
print(">>> m2 = Matrix('m2',[[1,2],[3,4]])")
print(">>> alpha = 1.0")
print(">>> m3 = et.dgemm(alpha,m1,m2)")
print(">>> print(m3)")
m1 = Matrix('m1',[[1,2],[3,4]])
m2 = Matrix('m2',[[1,2],[3,4]])
alpha = 1.0
m3 = et.dgemm(alpha,m1,m2)
print(m3)
print("\nScipy version")
print(">>> m1 = [[1,2],[3,4]]")
print(">>> m2 = [[1,2],[3,4]]")
print(">>> m3 = blas.dgemm(alpha,m1,m2)")
print(">>> print(m3)")
m1 = [[1,2],[3,4]]
m2 = [[1,2],[3,4]]
m3 = blas.dgemm(alpha,m1,m2)
print(m3)

print("\n(1b) - DGEMM")
print(">>> m1 = Matrix('m1',[[1,2,3],[4,5,6]])")
print(">>> m2 = Matrix('m2',[[1,2],[3,4],[5,6]])")
print(">>> alpha = 1.0")
print(">>> m3 = et.dgemm(alpha,m1,m2)")
print(">>> print(m3)")
m1 = Matrix('m1',[[1,2,3],[4,5,6]])
m2 = Matrix('m2',[[1,2],[3,4],[5,6]])
alpha = 1.0
m3 = et.dgemm(alpha,m1,m2)
print(m3)
print("\nScipy version")
print(">>> m1 = [[1,2,3],[4,5,6]]")
print(">>> m2 = [[1,2],[3,4],[5,6]]")
print(">>> m3 = blas.dgemm(alpha,m1,m2)")
print(">>> print(m3)")
m1 = [[1,2,3],[4,5,6]]
m2 = [[1,2],[3,4],[5,6]]
m3 = blas.dgemm(alpha,m1,m2)
print(m3)

print("\n(1b) - DGEMM")
print(">>> m1 = Matrix('m1',[[1,2,3],[4,5,6]])")
print(">>> m2 = Matrix('m2',[[1,2],[3,4],[5,6]])")
print(">>> alpha = 1.0")
print(">>> m3 = et.dgemm(alpha,m2,m1)")
print(">>> print(m3)")
m1 = Matrix('m1',[[1,2,3],[4,5,6]])
m2 = Matrix('m2',[[1,2],[3,4],[5,6]])
alpha = 1.0
m3 = et.dgemm(alpha,m2,m1)
print(m3)
print("\nScipy version")
print(">>> m1 = [[1,2,3],[4,5,6]]")
print(">>> m2 = [[1,2],[3,4],[5,6]]")
print(">>> m3 = blas.dgemm(alpha,m2,m1)")
print(">>> print(m3)")
m1 = [[1,2,3],[4,5,6]]
m2 = [[1,2],[3,4],[5,6]]
m3 = blas.dgemm(alpha,m2,m1)
print(m3)
