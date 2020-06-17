import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import cKDTree
from etraj.etraj import Grid
import time


g = Grid("10D",10)

N = 1000
#   generate random 2d data
x = np.random.normal(0,1,N)
y = np.random.normal(0,1,N)
data = np.vstack((x,y,x,y,y,x,y,x,x,y)).T

g.set_grid(data)
print("k neighbors search")
#   nearest neighbor test
start = time.time()
g.query_neighbors(10)
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
nano = [len(nano_neighbors_radius[i]) for i in range(len(nano_neighbors_radius))]
ckdt = [len(ckdtree_neighbors_radius[i]) for i in range(len(ckdtree_neighbors_radius))]
diff = [nano[i] - ckdt[i] for i in range(len(nano))]
print(sum(diff))


g = Grid("10D",10)
#   generate random 2d data

nano_times = []
ckdtree_times = []
nano_times_radius = []
ckdtree_times_radius = []
N = 10
nums = []
for i in range(10):
    x = np.random.normal(0,1,N)
    y = np.random.normal(0,1,N)
    data = np.vstack((x,y,x,y,x,y,x,y,x,y)).T
    g.set_grid(data)
    start = time.time()
    g.query_neighbors(10)
    end = time.time()
    nano_times.append(end - start)
    print("Nanoflann: ", end - start)
    nano_neighbors = g.get_neighbors()
    nano_distances = g.get_distances()
    start = time.time()
    tree = cKDTree(data)
    ckdtree_distances, ckdtree_neighbors = tree.query(data,10)
    end = time.time()
    ckdtree_times.append(end - start)
    print("cKDTree: ", end - start)

    print("\nRadius search")
    start = time.time()
    g.query_radius(2.0)
    end = time.time()
    nano_times_radius.append(end - start)
    print("Nanoflann: ", end - start)
    nano_neighbors_radius = g.get_neighbors_radius()
    nano_distances_radius = g.get_distances_radius()
    start = time.time()
    tree = cKDTree(data)
    ckdtree_neighbors_radius = tree.query_ball_point(data,2.0)
    end = time.time()
    ckdtree_times_radius.append(end - start)
    print("cKDTree: ", end - start)
    nums.append(N)
    N += 25
fig, axs = plt.subplots()
axs.plot(nums,nano_times,color='r',linestyle='--',label="NanoFLANN")
axs.plot(nums,ckdtree_times,color='k',linestyle='--',label="cKDTree")
axs.set_xlabel("N")
axs.set_ylabel("Time (sec)")
axs.set_title("Time (sec) vs. N, d = 10, k = 10")
plt.legend()
plt.savefig("s_vs_N_d10_k10.png")
plt.show()

fig, axs = plt.subplots()
axs.plot(nums,nano_times_radius,color='r',linestyle='--',label="NanoFLANN")
axs.plot(nums,ckdtree_times_radius,color='k',linestyle='--',label="cKDTree")
axs.set_xlabel("N")
axs.set_ylabel("Time (sec)")
axs.set_title("Time (sec) vs. N, d = 10, r = 2.0")
plt.savefig("s_vs_N_d10_r2.png")
plt.legend()
plt.show()
