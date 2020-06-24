import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix

#  tests for matrix class
print("\n-----------Tests for Matrix class--------------")
print("\n-----------------------------------------------")
print("\n(1) Constructors")

print("\n(1a) - Basic Constructor")
print(">>> m_1a = Matrix()")
print(">>> print(m_1a)")
m_1a = Matrix()
print(m_1a)

print("\n(1b) - (m x m) Constructor with m")
print(">>> m = 4")
print(">>> m_1b = Matrix(m)")
print(">>> print(m_1b)")
m = 4
m_1b = Matrix(m)
print(m_1b)

print("\n(1c) - (m x m) Constructor with m and name")
print(">>> m = 2")
print(">>> name = 'm_1c'")
print(">>> m_1c = Matrix(name,m)")
print(">>> print(m_1c)")
m = 2
name = 'm_1c'
m_1c = Matrix(name,m)
print(m_1c)

print("\n(1d) - (m x n) Constructor with m and n")
print(">>> m = 2")
print(">>> n = 3")
print(">>> m_1d = Matrix(m,n)")
print(">>> print(m_1d)")
m = 2
n = 3
m_1d = Matrix(m,n)
print(m_1d)

print("\n(1e) - (m x n) Constructor with name and m, n")
print(">>> m = 3")
print(">>> n = 2")
print(">>> name = 'm_1e'")
print(">>> m_1e = Matrix(name,m,n)")
print(">>> print(m_1e)")
m = 3
n = 2
name = 'm_1e'
m_1e = Matrix(name,m,n)
print(m_1e)

print("\n(1f) - (m x n) Constructor with m,n and initial value")
print(">>> m = 4")
print(">>> n = 2")
print(">>> init = 8.1")
print(">>> m_1f = Matrix(m,n,init)")
print(">>> print(m_1f)")
m = 4
n = 2
init = 8.1
m_1f = Matrix(m,n,init)
print(m_1f)

print("\n(1g) - (m x n) Constructor with name,m,n and initial value")
print(">>> name = 'm_1g'")
print(">>> m = 2")
print(">>> n = 6")
print(">>> init = 1.37")
print(">>> m_1g = Matrix(name,m,n,init)")
print(">>> print(m_1g)")
name = 'm_1g'
m = 2
n = 6
init = 1.37
m_1g = Matrix(name,m,n,init)
print(m_1g)

print("\n(1h) - (m x m) Constructor with m and flattened vector")
print(">>> m = 2")
print(">>> vec = [1,2,3,4]")
print(">>> m_1h = Matrix(m,vec)")
print(">>> print(m_1h)")
m = 2
vec = [1,2,3,4]
m_1h = Matrix(m,vec)
print(m_1h)

print("\n(1i) - (m x m) Constructor with name,m and flattened vector")
print(">>> name = 'm_1i'")
print(">>> m = 3")
print(">>> vec = [9,8,7,6,5,4,3,2,1]")
print(">>> m_1i = Matrix(name,m,vec)")
print(">>> print(m_1i)")
name = 'm_1i'
m = 3
vec = [9,8,7,6,5,4,3,2,1]
m_1i = Matrix(name,m,vec)
print(m_1i)

print("\n(1j) - (m x n) Constructor with m,n and flattened vector")
print(">>> m = 3")
print(">>> n = 2")
print(">>> vec = [1,2,3,4,5,6]")
print(">>> m_1j = Matrix(m,n,vec)")
print(">>> print(m_1j)")
m = 3
n = 2
vec = [1,2,3,4,5,6]
m_1j = Matrix(m,n,vec)
print(m_1j)

print("\n(1k) - (m x n) Constructor with name,m,n and flattened vector")
print(">>> name = 'm_1k'")
print(">>> m = 2")
print(">>> n = 3")
print(">>> vec = [8,7,6,5,4,3]")
print(">>> m_1k = Matrix(name,m,n,vec)")
print(">>> print(m_1k)")
name = 'm_1k'
m = 2
n = 3
vec = [8,7,6,5,4,3]
m_1k = Matrix(name,m,n,vec)
print(m_1k)

print("\n(1l) - (m x n) Constructor with 2d array")
print(">>> array = [[1,2],[3,4]]")
print(">>> m_1l = Matrix(array)")
print(">>> print(m_1l)")
array = [[1,2],[3,4]]
m_1l = Matrix(array)
print(m_1l)

print("\n(1m) - (m x n) Constructor with name and 2d array")
print(">>> name = 'm_1m'")
print(">>> array = [[1,2,3],[4,5,6]]")
print(">>> m_1m = Matrix(name,array)")
print(">>> print(m_1m)")
name = 'm_1m'
array = [[1,2,3],[4,5,6]]
m_1m = Matrix(name,array)
print(m_1m)

print("\n-----------------------------------------------")
print("\n(2) Getters and Setters")


print("\n-----------------------------------------------")
print("\n(3) Operator overloads")
