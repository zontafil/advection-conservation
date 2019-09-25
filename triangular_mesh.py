import numpy as np
import scipy
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from scipy.spatial import Voronoi, voronoi_plot_2d

# generate a triangular mesh with base b and nx, ny points
def generateMesh(b, nx, ny):
    points = np.array([[0, 0]])

    for i in range(nx):
        for j in range(ny):
            if (i==0) and (j == 0):
                continue
            offset = 0 if (i % 2) == 0 else (b/2)
            points = np.append(points, [[i*b, j*b + offset]], 0)

    return points