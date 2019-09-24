import numpy as np
import scipy
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d

b = 1.0
nx = 10
ny = 10

points = np.array([[0, 0]])

for i in range(nx):
    for j in range(ny):
        if (i==0) and (j == 0):
            continue
        offset = 0 if (i % 2) == 0 else (b/2)
        points = np.append(points, [[i, j + offset]], 0)

tri = Delaunay(points)
plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
plt.plot(points[:,0], points[:,1], 'o')
vor = Voronoi(points)
voronoi_plot_2d(vor)
plt.show()
