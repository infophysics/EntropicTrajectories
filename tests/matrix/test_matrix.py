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

print("\n(2a) - get_num_rows()")
print(">>> m_2a = Matrix(6,8)")
print(">>> m = m_2a.get_num_rows()")
print(">>> print(m)")
m_2a = Matrix(6,8)
m = m_2a.get_num_rows()
print(m)

print("\n(2b) - get_num_cols()")
print(">>> m_2b = Matrix(7,5)")
print(">>> n = m_2b.get_num_cols()")
print(">>> print(n)")
m_2b = Matrix(7,5)
n = m_2b.get_num_cols()
print(n)

print("\n(2c) - get_name()")
print(">>> m_2c = Matrix('m_2c',2)")
print(">>> name = m_2c.get_name()")
print(">>> print(name)")
m_2c = Matrix('m_2c',2)
name = m_2c.get_name()
print(name)

print("\n(2d) - get_array()")
print(">>> m_2d = Matrix([[1,2],[1,2]])")
print(">>> array = m_2d.get_array()")
print(">>> print(array)")
m_2d = Matrix([[1,2],[1,2]])
array = m_2d.get_array()
print(array)

print("\n(2e) - get_row(int k)")
print(">>> m_2e = Matrix([[1,2],[3,4]])")
print(">>> row1 = m_2e.get_row(0)")
print(">>> row2 = m_2d.get_row(1)")
print(">>> print(row1)")
print(">>> print(row2)")
m_2e = Matrix([[1,2],[3,4]])
row1 = m_2e.get_row(0)
row2 = m_2e.get_row(1)
print(row1)
print(row2)

print("\n(2f) - get_col(int k)")
print(">>> m_2f = Matrix([[1,2],[3,4]])")
print(">>> col1 = m_2f.get_col(0)")
print(">>> col2 = m_2d.get_col(1)")
print(">>> print(col1)")
print(">>> print(col2)")
m_2f = Matrix([[1,2],[3,4]])
col1 = m_2f.get_col(0)
col2 = m_2f.get_col(1)
print(col1)
print(col2)

print("\n(2g) - set_name(string name)")
print(">>> m_2g = Matrix(2,2)")
print(">>> name = 'm_2g'")
print(">>> m_2g.set_name(name)")
print(">>> print(m_2g)")
m_2g = Matrix(2,2)
name = 'm_2g'
m_2g.set_name(name)
print(m_2g)

print("\n(2h) - set_row(int k, [])")
print(">>> m_2h = Matrix(2,2)")
print(">>> row = [1,2]")
print(">>> m_2h.set_row(0,row)")
print(">>> print(m_2h)")
m_2h = Matrix(2,2)
row = [1,2]
m_2h.set_row(0,row)
print(m_2h)

print("\n(2i) - set_col(int k, [])")
print(">>> m_2i = Matrix(2,2)")
print(">>> col = [1,2]")
print(">>> m_2i.set_col(0,col)")
print(">>> print(m_2i)")
m_2i = Matrix(2,2)
col = [1,2]
m_2i.set_col(0,col)
print(m_2i)

print("\n(2j) - set_array(int m, [])")
print(">>> m_2j = Matrix()")
print(">>> array = [1,2,3,4]")
print(">>> m_2j.set_array(2,array)")
print(">>> print(m_2j)")
m_2j = Matrix()
array = [1,2,3,4]
m_2j.set_array(2,array)
print(m_2j)

print("\n(2k) - set_array([[]])")
print(">>> m_2k = Matrix()")
print(">>> array = [[1,2,3],[3,2,1]]")
print(">>> m_2k.set_array(array)")
print(">>> print(m_2k)")
m_2k = Matrix()
array = [[1,2,3],[3,2,1]]
m_2k.set_array(array)
print(m_2k)

print("\n-----------------------------------------------")
print("\n(3) Operator overloads")

print("\n(3a) Copy constructor")
print(">>> m1 = Matrix('m1',2,[0,4,2,1])")
print(">>> m2 = m1")
print(">>> print(m2)")
m1 = Matrix('m1',2,[0,4,2,1])
m2 = m1
print(m2)

print("\n(3b) Equivalence")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[1,2,3,4])")
print(">>> print(m1 =  m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[1,2,3,4])
print(m1 == m2)

print("\n(3c) Inequivalence")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[0,1,1,2])")
print(">>> print(m1 != m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[0,1,1,2])
print(m1 != m2)

print("\n(3d) Negation")
print(">>> m1 = Matrix('m1',2,[5,-6,7,0])")
print(">>> m2 = -m1")
print(">>> print(m2)")
m1 = Matrix('m1',2,[5,-6,7,0])
m2 = -m1
print(m2)

print("\n(3e) Matrix addition")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[4,3,2,1])")
print(">>> m3 = m1 + m2")
print(">>> print(m3)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[4,3,2,1])
m3 = m1 + m2
print(m3)

print("\n(3f) Matrix += addition")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[4,3,2,1])")
print(">>> m1 += m2")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[4,3,2,1])
m1 += m2
print(m1)

print("\n(3g) Matrix subtraction")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[4,3,2,1])")
print(">>> m3 = m1 - m2")
print(">>> print(m3)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[4,3,2,1])
m3 = m1 - m2
print(m3)

print("\n(3h) Matrix += subtraction")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[4,3,2,1])")
print(">>> m1 -= m2")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[4,3,2,1])
m1 -= m2
print(m1)

print("\n(3i) Scalar left-addition")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = m1 + 1.0")
print(">>> print(m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = m1 + 1.0
print(m2)

print("\n(3j) Scalar right-addition")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = 1.0 + m1")
print(">>> print(m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = 1.0 + m1
print(m2)

print("\n(3k) Scalar left-subtraction")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = m1 - 1.0")
print(">>> print(m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = m1 - 1.0
print(m2)

print("\n(3k) Scalar right-subtraction")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = 1.0 - m1")
print(">>> print(m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = 1.0 - m1
print(m2)

print("\n(3l) Scalar left-multiplication")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = 3.14 * m1")
print(">>> print(m2)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = 3.14 * m1
print(m2)

print("\n(3m) Scalar right-multiplication")
print(">>> m1 = Matrix('m1',2,3,[1,2,3,4,5,6])")
print(">>> m2 = m1 * 3.14")
print(">>> print(m2)")
m1 = Matrix('m1',2,3,[1,2,3,4,5,6])
m2 = m1 * 3.14
print(m2)

print("\n(3n) Scalar left-division")
print(">>> m1 = Matrix('m1',3,2,[1,2,3,4,5,6])")
print(">>> m2 = 3.14 / m1")
print(">>> print(m2)")
m1 = Matrix('m1',3,2,[1,2,3,4,5,6])
m2 = 3.14 / m1
print(m2)

print("\n(3o) Scalar right-division")
print(">>> m1 = Matrix('m1',2,3,[1,2,3,4,5,6])")
print(">>> m2 = m1 / 3.14")
print(">>> print(m2)")
m1 = Matrix('m1',2,3,[1,2,3,4,5,6])
m2 = m1 / 3.14
print(m2)

print("\n(3p) Scalar inplace addition")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m1 += 3.145")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m1 += 3.145
print(m1)

print("\n(3q) Scalar inplace subtraction")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m1 -= 3.145")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m1 -= 3.145
print(m1)

print("\n(3r) Scalar inplace multiplication")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m1 *= 3.145")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m1 *= 3.145
print(m1)

print("\n(3s) Scalar inplace division")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m1 /= 3.145")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m1 /= 3.145
print(m1)

print("\n(3t) Const access")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> x = m1[0,1]")
print(">>> print(x)")
m1 = Matrix('m1',2,[1,2,3,4])
x = m1[0,1]
print(x)

print("\n(3u) Assignment")
print(">>> m1 = Matrix('m1',2,3,[1,2,3,4,5,6])")
print(">>> m1[0,1] = 7.8")
print(">>> print(m1)")
m1 = Matrix('m1',2,3,[1,2,3,4,5,6])
m1[0,1] = 7.8
print(m1)

print("\n(3v) Multiplication")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[4,3,2,1])")
print(">>> m3 = m1 * m2")
print(">>> print(m3)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[4,3,2,1])
m3 = m1 * m2
print(m3)

print("\n(3w) Inplace Multiplication")
print(">>> m1 = Matrix('m1',2,[1,2,3,4])")
print(">>> m2 = Matrix('m2',2,[4,3,2,1])")
print(">>> m1 *= m2")
print(">>> print(m1)")
m1 = Matrix('m1',2,[1,2,3,4])
m2 = Matrix('m2',2,[4,3,2,1])
m1 *= m2
print(m1)

print("\n-----------------------------------------------")
print("\n(4) Various methods")
