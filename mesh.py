import numpy as np
import math
from scipy.spatial import Voronoi, Delaunay
from mesh_utils import inside_convex_polygon
from typing import List


class Mesh:
    """Delaunay and Voronoi meshes struct

    Attributes:
        delaunay {scipy.spatial.Delaunay} -- delaunay mesh
        voronoi {scipy.spatial.Voronoi} -- voronoi mesh
        voronoiNeighbors {list} -- list of neighbors for each voronoi region
        voronoiRegions {list} -- list of voronoiregions, cleand from the python library
    """
    def __init__(self, points: np.array):
        """

        Arguments:
            points {np.array} -- points which define the mesh
        """
        self.delaunay = Delaunay(points)
        self.voronoi = Voronoi(points)

        self.buildVoronoiRegions()
        self.buildVoronoiNeighborsIndex()

    def isPointInsideVoronoiRegion(self, point: List[float],
                                   regionIdx: int):
        """Check if a point is inside a given voronoi region

        Arguments:
            point {list} -- coordinates of the point
            regionIdx {int} -- index of the voronoi region

        Returns:
            bool -- True if the point is inside the region
        """
        # find region's vertices coordinates
        coordinates = []
        for i in self.voronoiRegions[regionIdx]:
            coordinates.append(self.voronoi.vertices[i])
        return inside_convex_polygon(point, coordinates)

    def buildVoronoiRegions(self):
        """Clean the voronoi regions array and get rid of void elements
        """
        self.voronoiRegions = []
        for i in self.voronoi.regions:
            if (len(i) > 0):
                self.voronoiRegions.append(i)

    def buildVoronoiNeighborsIndex(self):
        """Build An array of neighbors for each voronoi region
           populate the self.voronoiNeighbors array
        """
        self.voronoiNeighbors = [None]*len(self.voronoiRegions)
        for i, region_i in enumerate(self.voronoiRegions):
            self.voronoiNeighbors[i] = []
            for j, region_j in enumerate(self.voronoiRegions):
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


class EquilateralTriangularMesh(Mesh):
    def __init__(self, x0, x1, y0, y1, b):
        """Generate a equilateral triangular mesh and its dual

        Arguments:
            x0 {float} -- lower x limit of the mesh
            x1 {float} -- upper x limit of the mesh
            y0 {float} -- lower y limit of the mesh
            y1 {float} -- upper y limit of the mesh
            b {float} -- base value of the triangles

        Returns:
            Mesh -- return the primary (Delaunay) and dual (Voronoi) meshes
        """
        points = np.array([[]])

        npAxis = 1
        i = 0
        px = x0
        py = y0
        while (py < y1):
            while (px < x1):
                points = np.append(points, [[px, py]], npAxis)
                npAxis = 0
                px += b
            i = i + 1
            px = x0 if (i % 2) == 0 else x0 + (b/2)
            py = py + b * math.sqrt(3) / 2.0

        super().__init__(points)
