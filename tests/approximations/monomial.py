import numpy as np
import matplotlib.pyplot as plt
from etraj.etraj import Matrix
import etraj.etraj as et

xyz = et.monomial_n(3,3)
print(xyz)

x = et.taylor_monomial_factors(3,3)
print(x)

z = et.taylor_polynomial(1,1.5,3)
print(z)

y = et.taylor_monomial_expansion([1,1,1],[1.5,1.5,1.5],3)
print(y)
