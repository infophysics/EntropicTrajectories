import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix, UGrid, Approximator, ScalarField
import etraj.etraj as et
from scipy import linalg
from scipy.linalg import svdvals

#  tests for matrix class
print("\n--------Tests for Approximator class-----------")
print("\n-----------------------------------------------")
print("\n(1) Constructors")

print("\n(1a) - Basic constructor")
print(">>> app = Approximator()")
print(">>> print(app)")
app = Approximator()
print(app)

print("\n(1b) - Constructor with int type")
print(">>> app = Approximator(1)")
print(">>> print(app)")
app = Approximator(1)
print(app)

print("\n(1c) - Constructor with string type")
print(">>> app = Approximator('MLS')")
print(">>> print(app)")
app = Approximator('MLS')
print(app)

print("\n(2) Getters and setters")

print("\n(2a) - ApproxType")
print(">>> app = Approximator()")
print(">>> type = app.get_approx_type()")
print(">>> print(type)")
app = Approximator()
type = app.get_approx_type()
print(type)

print("\n(2b) - ApproxParams")
print(">>> app = Approximator()")
print(">>> app.get_approx_params().k = 5")
print(">>> print(app)")
app = Approximator()
a = app.get_approx_params()
a.k = 4
app.set_approx_params(a)
print(app)
