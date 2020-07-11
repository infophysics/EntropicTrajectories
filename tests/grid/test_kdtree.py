import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from etraj.etraj import UGrid
import time




N = 1000
#   generate random 2d data
x = np.random.normal(0,1,N)
y = np.random.normal(0,1,N)
data = np.vstack((x,y)).T

g = UGrid(data)

g.query_neighbors(5)
neighbors = g.get_neighbors()
#neighbors2 = g.query_neighbors(x,5)
neighbor = g.query_neighbors([0.46,0.46],3)
