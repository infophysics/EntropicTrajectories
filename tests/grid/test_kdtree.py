import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from etraj.etraj import UGrid
import time


g = UGrid("10D",10)

N = 1000
#   generate random 2d data
x = np.random.normal(0,1,N)
y = np.random.normal(0,1,N)
data = np.vstack((x,y,x,y,y,x,y,x,x,y)).T

g.set_ugrid(data)

g.query_neighbors(5)
neighbors = g.get_neighbors()

neighbors2 = g.query_neighbors(data,5)

vals = [[neighbors[i][j] - neighbors2[i][j] for j in range(len(neighbors[i]))] for i in range(len(neighbors))]
print(vals)
