import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from etraj.etraj import Grid
import time


g = Grid("10D",10)

N = 10
#   generate random 2d data
x = np.random.normal(0,1,N)
y = np.random.normal(0,1,N)
data = np.vstack((x,y,x,y,y,x,y,x,x,y)).T

g.set_grid(data)
print("k neighbors search")
#   nearest neighbor test
start = time.time()
g.query(10)
end = time.time()
print("Nanoflann: ", end - start)
nano_neighbors = g.get_neighbors()
nano_distances = g.get_distances()
start = time.time()
tree = cKDTree(data)
ckdtree_distances, ckdtree_neighbors = tree.query(data,10)
end = time.time()
print("cKDTree: ", end - start)
print(nano_neighbors == ckdtree_neighbors)
print(np.sqrt(nano_distances) - ckdtree_distances)

print("\nRadius search")
start = time.time()
g.query_radius(2.0)
end = time.time()
print("Nanoflann: ", end - start)
nano_neighbors_radius = g.get_neighbors_radius()
nano_distances_radius = g.get_distances_radius()
start = time.time()
tree = cKDTree(data)
ckdtree_neighbors_radius = tree.query_ball_point(data,2.0)
end = time.time()
print("cKDTree: ", end - start)
print(nano_neighbors_radius)
print(ckdtree_neighbors_radius)
