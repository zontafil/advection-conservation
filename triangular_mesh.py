import numpy as np
from scipy.spatial import Voronoi, Delaunay


class Mesh:
    """Delaunay and Voronoi meshes struct

    Attributes:
        delaunay {scipy.spatial.Delaunay} -- delaunay mesh
        voronoi {scipy.spatial.Voronoi} -- voronoi mesh
        voronoiNeighbors {list} -- list of neighbors for each voronoi region
    """
    def __init__(self, delaunay: Delaunay, voronoi: Voronoi):
        """

        Arguments:
            delaunay {scipy.spatial.Delaunay} -- delaunay mesh
            voronoi {scipy.spatial.Voronoi} -- voronoi mesh
        """
        self.delaunay = delaunay
        self.voronoi = voronoi

    def buildVoronoiNeighborsIndex(self):
        """Build An array of neighbors for each voronoi region
           populate the self.voronoiNeighbors array
        """
        self.voronoiNeighbors = [None]*len(self.voronoi.regions)
        for i, region_i in enumerate(self.voronoi.regions):
            self.voronoiNeighbors[i] = []
            for j, region_j in enumerate(self.voronoi.regions):
                if (i == j):
                    continue
                # check if region i and j are neighbors.
                # Add them to the list in case
                if self.checkRegionsHaveCommonEdge(region_i, region_j):
                    if (j not in self.voronoiNeighbors[i]):
                        self.voronoiNeighbors[i].append(j)

    def findVoronoiRegionsNeighbors(self, regionsIndexes):
        """Find all the neighbors of a set of voronoi mesh elements (a.k.a regions)

        Arguments:
            regionsIndexes {list} -- indexes (int) of the region
        """
        neighbors = []
        for i in regionsIndexes:
            for neighbor_i in self.voronoiNeighbors[i]:
                if (neighbor_i not in neighbors):
                    neighbors.append(neighbor_i)
        return neighbors

    def checkRegionsHaveCommonEdge(self, a, b):
        """Check if two regions have at least one common edge

        Arguments:
            a {list} -- region A object
            b {list} -- region B object

        Returns:
            True/False
        """
        ret = False
        for i in a:
            if i == -1:
                continue
            for j in b:
                if j == -1:
                    continue
                if i == j:
                    ret = True
        return ret


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

    mesh = Mesh(Delaunay(points), Voronoi(points))

    # precompute useful stuff, i.e. voronoi neighbors index
    mesh.buildVoronoiNeighborsIndex()
    return mesh
