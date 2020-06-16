import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix

#  tests for matrix class
print("...Matrix tests...")

print("(1) Constructors...")
print("(1a) - m_1a = Matrix(4)")
m_1a = Matrix(4)
print(m_1a)

print("(1b) - m_1b = Matrix('m_1b', 4)")
m_1b = Matrix('m_1b',4)
print(m_1b)

print("(1c) - m_1c = Matrix(2, 3)")
m_1c = Matrix(2,3)
print(m_1c)

print("(1d) - m_1d = Matrix('m_1d', 3, 2)")
m_1d = Matrix('m_1d',3,2)
print(m_1d)

print("(1e) - m_1e = Matrix(2, 3, 3.14)")
m_1e = Matrix(2,3,3.14)
print(m_1e)

print("(1f) - m_1f = Matrix('m_1f', 3, 2, 3.14)")
m_1f = Matrix('m_1f',3,2,3.14)
print(m_1f)

print("(1g) - m_1g = Matrix(2, [4,3,2,1])")
m_1g = Matrix(2, [4,3,2,1])
print(m_1g)

print("(1h) - m_1h = Matrix('m_1h', 2, [1,2,3,4])")
m_1h = Matrix('m_1h', 2, [1,2,3,4])
print(m_1h)

print("(1i) - m_1i = Matrix(2, 4, [1,2,3,4,5,6,7,8])")
m_1i = Matrix(2, 4, [1,2,3,4,5,6,7,8])
print(m_1i)

print("(1j) - m_1j = Matrix('m_1j', 1, 3, [2,4,6])")
m_1j = Matrix('m_1j', 1, 3, [2,4,6])
print(m_1j)

print("(1k) - m_1k = Matrix([[1,3],[5,6]])")
m_1k = Matrix([[1,3],[5,6]])
print(m_1k)

print("(1l) - m_1l = Matrix('m_1l', [[1,3],[5,6]])")
m_1l = Matrix('m_1l', [[1,3],[5,6]])
print(m_1l)

print("(2) Getters and setters...")
print("(2a) - x = m_1l.get_rows()")
x = m_1l.get_rows()
print("Number of rows in m_1l: ", x)

print("(2b) - x = m_1l.get_cols()")
x = m_1l.get_cols()
print("Number of cols in m_1l: ", x)

print("(2c) - n = m_1l.get_name()")
n = m_1l.get_name()
print("Name of m_1l: ", n)

print("(2d) - mat = m_1l.get_mat()")
mat = m_1l.get_mat()
print("Flattened matrix of m_1l: ", mat)

print("(2e) - row = m_1l.get_row(0)")
row = m_1l.get_row(0)
print("Row 0 of m_1l: ", row)

print("(2f) - col = m_1l.get_col(0)")
col = m_1l.get_col(0)
print("Column 0 of m_1l: ", col)


print("(3) Operator overloads...")


m2 = Matrix("B",2,2,[1,2,3,4])

m3 = m_1a * m2
m3.set_name("C")
print(m3)

m4 = Matrix("D",3,10,[np.random.normal(0,10,1)[0] for i in range(100)])

print(m4)
m5 = Matrix("E",10,3,[np.random.normal(0,10,1)[0] for i in range(100)])
print(m5)

m6 = Matrix(np.asarray([[1,2],[1,2]]))
print(m6)

print(m6)

print(m6)
m7 = 0.0 / m6
print(m7)

m7[0,0] = 8.0
