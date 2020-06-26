import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Vector

#  tests for vector class
print("\n-----------Tests for Vector class--------------")
print("\n-----------------------------------------------")
print("\n(1) Constructors")

print("\n(1a) Basic Constructor")
print(">>> v = Vector()")
print(">>> print(v)")
v = Vector()
print(v)

print("\n(1b) Constructor with dim")
print(">>> dim = 6")
print(">>> v = Vector(dim)")
print(">>> print(v)")
dim = 6
v = Vector(dim)
print(v)

print("\n(1c) Constructor with name and dim")
print(">>> name = 'v'")
print(">>> dim = 5")
print(">>> v = Vector(name,dim)")
print(">>> print(v)")
name = 'v'
dim = 5
v = Vector(name,dim)
print(v)

print("\n(1d) Constructor with vector")
print(">>> vec = [0,1,2,3,4,5,6]")
print(">>> v = Vector(vec)")
print(">>> print(v)")
vec = [0,1,2,3,4,5,6]
v = Vector(vec)
print(v)

print("\n(1e) Constructor with vec and name")
print(">>> name = 'vec_with_name'")
print(">>> vec = [6,5,4,3,2,1]")
print(">>> v = Vector(name,vec)")
print(">>> print(v)")
name = 'vec_with_name'
vec = [6,5,4,3,2,1]
v = Vector(name,vec)
print(v)

print("\n(1f) Constructor with dim and init value")
print(">>> dim = 10")
print(">>> init = 3.14567")
print(">>> v = Vector(dim,init)")
print(">>> print(v)")
dim = 10
init = 3.141567
v = Vector(dim,init)
print(v)

print("\n(1g) Constructor with name, dim and init value")
print(">>> name = 'v_with_name_and_init'")
print(">>> dim = 3")
print(">>> init = 1.27e-5")
print(">>> v = Vector(name,dim,init)")
print(">>> print(v)")
name = 'v_with_name_and_init'
dim = 3
init = 1.27e-5
v = Vector(name,dim,init)
print(v)

print("\n-----------------------------------------------")
print("\n(2) Getters and Setters")

print("\n(2a) Get dim")
print(">>> v = Vector('v',[4,3,2,6])")
print(">>> dim = v.get_dim()")
print(">>> print(dim)")
v = Vector('v',[4,3,2,6])
dim = v.get_dim()
print(dim)

print("\n(2b) Get vector")
print(">>> v = Vector('random',np.random.normal(0,1,14))")
print(">>> x = v.get_vec()")
print(">>> print(x)")
v = Vector('random',np.random.normal(0,1,14))
x = v.get_vec()
print(x)

print("\n(2c) Get name")
print(">>> v = Vector('special_name',6)")
print(">>> name = v.get_name()")
print(">>> print(name)")
v = Vector('special_name',6)
name = v.get_name()
print(name)

print("\n(2d) Set dim")
print(">>> v = Vector()")
print(">>> v.set_dim(10)")
print(">>> print(v)")
v = Vector()
v.set_dim(10)
print(v)

print("\n(2e) Set vector")
print(">>> v = Vector('vec',5)")
print(">>> x = [4,5,3]")
print(">>> v.set_vec(x)")
print(">>> print(v)")
v = Vector('vec',5)
x = [4,5,3]
v.set_vec(x)
print(v)

print("\n(2f) Set name")
print(">>> v = Vector('orig_name',5)")
print(">>> v.set_name('new_name')")
print(">>> print(v)")
v = Vector('orig_name',5)
v.set_name('new_name')
print(v)

print("\n-----------------------------------------------")
print("\n(3) Operator overloads")

print("\n(3a) Copy constructor")
print(">>> v = Vector('v',[0,4,2,1])")
print(">>> u = v")
print(">>> print(u)")
v = Vector('v',[0,4,2,1])
u = v
print(u)

print("\n(3b) Equivalence")
print(">>> v = Vector('v',[1,2,3,4])")
print(">>> u = Vector('u',[1,2,3,4])")
print(">>> print(v == u)")
v = Vector('v',[1,2,3,4])
u = Vector('u',[1,2,3,4])
print(v == u)

print("\n(3c) Inequivalence")
print(">>> v = Vector('v',[1,2,3])")
print(">>> u = Vector('u',[0,1])")
print(">>> print(v != u)")
v = Vector('v',[1,2,3])
u = Vector('u',[0,1])
print(v != 0)

print("\n(3d) Negation")
print(">>> v = Vector('v',[5,-6,7,0])")
print(">>> u = -v")
print(">>> print(u)")
v = Vector('v',[5,-6,7,0])
u = -v
print(u)

print("\n(3e) Vector addition")
print(">>> v = Vector('v',[1,2,3])")
print(">>> u = Vector('u',[3,2,1])")
print(">>> w = v + u")
print(">>> print(w)")
v = Vector('v',[1,2,3])
u = Vector('u',[3,2,1])
w = v + u
print(w)

print("\n(3f) Vector += addition")
print(">>> v = Vector('v',[1,2,3])")
print(">>> u = Vector('u',[3,2,1])")
print(">>> v += u")
print(">>> print(v)")
v = Vector('v',[1,2,3])
u = Vector('u',[3,2,1])
v += u
print(v)

print("\n(3g) Vector subtraction")
print(">>> v = Vector('v',[1,2,3])")
print(">>> u = Vector('u',[3,2,1])")
print(">>> w = v - u")
print(">>> print(w)")
v = Vector('v',[1,2,3])
u = Vector('u',[3,2,1])
w = v - u
print(w)

print("\n(3h) Vector += subtraction")
print(">>> v = Vector('v',[1,2,3])")
print(">>> u = Vector('u',[3,2,1])")
print(">>> v -= u")
print(">>> print(v)")
v = Vector('v',[1,2,3])
u = Vector('u',[3,2,1])
v -= u
print(v)

print("\n(3i) Scalar left-addition")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> u = v + 1.0")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5])
u = v + 1.0
print(u)

print("\n(3j) Scalar right-addition")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> u = 1.0 + v")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5])
u = 1.0 + v
print(u)

print("\n(3k) Scalar left-subtraction")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> u = v - 1.0")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5])
u = v - 1.0
print(u)

print("\n(3k) Scalar right-subtraction")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> u = 1.0 - v")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5])
u = 1.0 - v
print(u)

print("\n(3l) Scalar left-multiplication")
print(">>> v = Vector('v',[1,2,3,4,5,6])")
print(">>> u = 3.14 * v")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5,6])
u = 3.14 * v
print(u)

print("\n(3m) Scalar right-multiplication")
print(">>> v = Vector('v',[1,2,3,4,5,6])")
print(">>> u = v * 3.14")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5,6])
u = v * 3.14
print(u)

print("\n(3n) Scalar left-division")
print(">>> v = Vector('v',[1,2,3,4,5,6])")
print(">>> u = 3.14 / v")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5,6])
u = 3.14 / v
print(u)

print("\n(3o) Scalar right-division")
print(">>> v = Vector('v',[1,2,3,4,5,6])")
print(">>> u = v / 3.14")
print(">>> print(u)")
v = Vector('v',[1,2,3,4,5,6])
u = v / 3.14
print(u)

print("\n(3p) Scalar inplace addition")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> v += 3.145")
print(">>> print(v)")
v = Vector('v',[1,2,3,4,5])
v += 3.145
print(v)

print("\n(3q) Scalar inplace subtraction")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> v -= 3.145")
print(">>> print(v)")
v = Vector('v',[1,2,3,4,5])
v -= 3.145
print(v)

print("\n(3r) Scalar inplace multiplication")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> v *= 3.145")
print(">>> print(v)")
v = Vector('v',[1,2,3,4,5])
v *= 3.145
print(v)

print("\n(3s) Scalar inplace division")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> v /= 3.145")
print(">>> print(v)")
v = Vector('v',[1,2,3,4,5])
v /= 3.145
print(v)

print("\n(3t) Const access")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> x = v[3]")
print(">>> print(x)")
v = Vector('v',[1,2,3,4,5])
x = v[3]
print(x)

print("\n(3u) Assignment")
print(">>> v = Vector('v',[1,2,3,4,5])")
print(">>> v[3] = 7.8")
print(">>> print(v)")
v = Vector('v',[1,2,3,4,5])
v[3] = 7.8
print(v)

print("\n-----------------------------------------------")
print("\n(4) Various geometric methods")

print("\n(4a) Dot product")
print(">>> v = Vector('v',[1,2,3,4])")
print(">>> u = Vector('u',[4,3,2,1])")
print(">>> s = v.dot(u)")
print(">>> print(s)")
v = Vector('v',[1,2,3,4])
u = Vector('u',[4,3,2,1])
s = v.dot(u)
print(s)

print("\n(4b) Dot product operator overload")
print(">>> v = Vector('v',[1,2,3,4])")
print(">>> u = Vector('u',[4,3,2,1])")
print(">>> s = v*u")
print(">>> print(s)")
v = Vector('v',[1,2,3,4])
u = Vector('u',[4,3,2,1])
s = v*u
print(s)
