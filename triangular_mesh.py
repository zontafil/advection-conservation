import numpy as np
import scipy


class Mesh:
    """Delaunay and Voronoi meshes struct

    Attributes:
        delaunay {scipy.spatial.Delaunay} -- delaunay mesh
        voronoi {scipy.spatial.Voronoi} -- voronoi mesh
    """
    def __init__(self, delaunay, voronoi):
        """

        Arguments:
            delaunay {scipy.spatial.Delaunay} -- delaunay mesh
            voronoi {scipy.spatial.Voronoi} -- voronoi mesh
        """
        self.delaunay = delaunay
        self.voronoi = voronoi


def generateMesh(b, nx, ny):
    """Generate a equilateral triangular mesh and its dual

    Arguments:
        b {int} -- base length of each triangle
        nx {int} -- number of points in the x dimension
        ny {int} -- number of points in the y dimension

    Returns:
        Mesh -- return the primary (Delaunay) and dual (Voronoi) meshes
    """
    points = np.array([[0, 0]])

    for i in range(nx):
        for j in range(ny):
            # FIXME - remove this ugly check
            if (i == 0) and (j == 0):
                continue
            offset = 0 if (i % 2) == 0 else (b/2)
            points = np.append(points, [[i*b, j*b + offset]], 0)

    return Mesh(scipy.spatial.Delaunay(points), scipy.spatial.Voronoi(points))
