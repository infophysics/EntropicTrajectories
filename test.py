import numpy as np
import matplotlib.pyplot as plt

from etraj.etraj import Matrix

m1 = Matrix("A",2,2,[3,1,-1,2])
m2 = Matrix("B",2,2,[1,2,3,4])

m3 = m1 * m2
m3.set_name("C")
print(m3)
