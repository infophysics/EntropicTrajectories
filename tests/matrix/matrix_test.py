import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix

#  tests for matrix class
print("...Matrix tests...")

print("\n(1) Constructors...")
print("\n(1a) - m_1a = Matrix(4)")
m_1a = Matrix(4)
print(m_1a)

print("\n(1b) - m_1b = Matrix('m_1b', 4)")
m_1b = Matrix('m_1b',4)
print(m_1b)

print("\n(1c) - m_1c = Matrix(2, 3)")
m_1c = Matrix(2,3)
m_1c.set_name('m_1c')
print(m_1c)

print("\n(1d) - m_1d = Matrix('m_1d', 3, 2)")
m_1d = Matrix('m_1d',3,2)
print(m_1d)

print("\n(1e) - m_1e = Matrix(2, 3, 3.14)")
m_1e = Matrix(2,3,3.14)
print(m_1e)

print("\n(1f) - m_1f = Matrix('m_1f', 3, 2, 3.14)")
m_1f = Matrix('m_1f',3,2,3.14)
print(m_1f)

print("\n(1g) - m_1g = Matrix(2, [4,3,2,1])")
m_1g = Matrix(2, [4,3,2,1])
m_1g.set_name('m_1g')
print(m_1g)

print("\n(1h) - m_1h = Matrix('m_1h', 2, [1,2,3,4])")
m_1h = Matrix('m_1h', 2, [1,2,3,4])
print(m_1h)

print("\n(1i) - m_1i = Matrix(2, 4, [1,2,3,4,5,6,7,8])")
m_1i = Matrix(2, 4, [1,2,3,4,5,6,7,8])
print(m_1i)

print("\n(1j) - m_1j = Matrix('m_1j', 1, 3, [2,4,6])")
m_1j = Matrix('m_1j', 1, 3, [2,4,6])
print(m_1j)

print("\n(1k) - m_1k = Matrix([[1,3],[5,6]])")
m_1k = Matrix([[1,3],[5,6]])
print(m_1k)

print("\n(1l) - m_1l = Matrix('m_1l', [[1,3],[5,6]])")
m_1l = Matrix('m_1l', [[1,3],[5,6]])
print(m_1l)

print("\n(2) Getters and setters...")
print("\n(2a) - x = m_1l.get_num_rows()")
x = m_1l.get_num_rows()
print("Number of rows in m_1l: ", x)

print("\n(2b) - x = m_1l.get_num_cols()")
x = m_1l.get_num_cols()
print("Number of cols in m_1l: ", x)

print("\n(2c) - n = m_1l.get_name()")
n = m_1l.get_name()
print("Name of m_1l: ", n)

print("\n(2d) - mat = m_1l.get_array()")
mat = m_1l.get_array()
print("Flattened matrix of m_1l: ", mat)

print("\n(2e) - row = m_1l.get_row(0)")
row = m_1l.get_row(0)
print("Row 0 of m_1l: ", row)

print("\n(2f) - col = m_1l.get_col(0)")
col = m_1l.get_col(0)
print("Column 0 of m_1l: ", col)


print("\n(3) Operator overloads...")
print("\n(3a) - m_1a == m_2a")
print(m_1a == m_1b)

print("\n(3a) - m_1g + m_2h")
print(m_1g + m_1h)

print("\n(3a') - m_1h + m_2g")
print(m_1h + m_1g)

print("\n(3b) - m_1g - m_2h")
print(m_1g - m_1h)

print("\n(3b') - m_1h - m_2g")
print(m_1h - m_1g)

print("\n(3c) - m_1h += m_2g")
m_1h += m_1g
print(m_1h)

print("\n(3d) - m_1g -= m_2h")
m_1g -= m_1h
print(m_1g)

print("\n(3e) - m_3a * m_3b")
m_3a = Matrix('m_3a', 2, [4,3,2,1])
m_3b = Matrix('m_3b', 2, [1,2,3,4])
m_3ap = m_3a * m_3b
m_3bp = m_3b * m_3a
print(m_3ap)
print(m_3bp)

print("\n(3e) - m_3b * m_3a")
print(m_3bp)

print("\n(3f) - m_3a *= m_3b")
m_3a *= m_3b
print(m_3a)

print("\n(3g) - m_3a + 3.14")
m_3a = Matrix('m_3a', 2, [4,3,2,1])
m_3b = m_3a + 3.14
print(m_3b)

print("\n(3g') - 3.14 + m_3a")
m_3b = 3.14 + m_3a
print(m_3b)

print("\n(3h) - m_3a - 3.14")
m_3a = Matrix('m_3a', 2, [4,3,2,1])
m_3b = m_3a - 3.14
print(m_3b)

print("\n(3h') - 3.14 - m_3a")
m_3b = 3.14 - m_3a
print(m_3b)

print("\n(3i) - m_3a * 3.14")
m_3a = Matrix('m_3a', 2, [4,3,2,1])
m_3b = m_3a * 3.14
print(m_3b)

print("\n(3i') - 3.14 * m_3a")
m_3b = 3.14 * m_3a
print(m_3b)

print("\n(3j) - m_3a / 3.14")
m_3a = Matrix('m_3a', 2, [4,3,2,1])
m_3b = m_3a / 3.14
print(m_3b)

print("\n(3j') - 3.14 / m_1a")
m_3b = 3.14 / m_3a
print(m_3b)

print("\n(4) Matrix access operator...")

print("\n(4a) - m_4a[0,0] = 3.14")
m_4a = Matrix('m_4a', 2, [4,3,2,1])
m_4a[0,0] = 3.14
print(m_4a)



x = np.asarray(m_4a)
print(x)
