import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix, UGrid, Interpolator, ScalarField
import etraj.etraj as et
from scipy import linalg
from scipy.linalg import svdvals

#  tests for matrix class
print("\n--------Tests for Interpolator class-----------")
print("\n-----------------------------------------------")
print("\n(1) Constructors")

print("\n(1a) - Basic constructor")
print(">>> app = Interpolator()")
print(">>> print(app)")
app = Interpolator()
print(app)

print("\n(1b) - Constructor with int type")
print(">>> app = Interpolator(1)")
print(">>> print(app)")
app = Interpolator(1)
print(app)

print("\n(1c) - Constructor with string type")
print(">>> app = Interpolator('MLS')")
print(">>> print(app)")
app = Interpolator('MLS')
print(app)

print("\n(2) Getters and setters")

print("\n(2a) - InterpolatorType")
print(">>> app = Interpolator()")
print(">>> type = app.get_approx_type()")
print(">>> print(type)")
app = Interpolator()
type = app.get_approx_type()
print(type)

print("\n(2b) - InterpolatorParams")
print(">>> app = Interpolator()")
print(">>> app.get_approx_params().k = 5")
print(">>> print(app)")
app = Interpolator()
a = app.get_approx_params()
a.k = 4
app.set_approx_params(a)
print(app)
