from mesh import Mesh
from singleParticleIntegrator import Particle
import numpy as np


def computeBarycentric(p: list, Q: list):
    """Compute the barycentric weights for a point p in an n-gon Q
        Assumes p is strictly within Q and the vertices qj of Q are ordered

    Arguments:
        p {list} -- coordinates of the point inside the polygon
        Q {list} -- coordinates of the polygon vertices
    Returns:
        [list] -- weights for each polygon's vertex
    """
    n = len(Q)
    w = [None]*len(Q)
    weightSum = 0
    for j, qj in enumerate(Q):
        prev = (j + n - 1) % n
        next = (j + 1) % n
        w[j] = (cotangent(p, qj, Q[prev]) + cotangent(p, qj, Q[next])) / np.linalg.norm(np.subtract(p, Q[j])**2)
        weightSum += w[j]

    # normalize by the area of the polygon
    return np.divide(w, weightSum)


def cotangent(a: list, b: list, c: list):
    """Compute the cotangent of the non-degenerate triangle abc at vertex b

    Arguments:
        a {list} -- coordinates of the triangle's point a
        b {list} -- coordinates of the triangle's point b
        c {list} -- coordinates of the triangle's point c

    Returns:
        [float]
    """
    ba = np.subtract(a, b)
    bc = np.subtract(c, b)
    return (np.dot(bc, ba) / abs(np.cross(bc, ba)))


class ParticlesMeshDeposit:
    """Deposit and distribution function matrices and methods

    Attributes:
        mesh
        f -- distribution function

    """
    def __init__(self, mesh: Mesh):
        """Constructor. A valid mesh must be passed

        Arguments:
            mesh {Mesh}
        """
        self.mesh = mesh
        self.f: list(float) = [0]*len(mesh.delaunay.simplices)

    def addParticle(self, particle: Particle, lastKnownVoronoiPosition=0):
        """Add a particle and deposit it to the dual mesh and to the distribution function

        Arguments:
            particle {Particle}

        Keyword Arguments:
            lastKnownVoronoiPosition {int} -- Last known voronoi region position of the particle.
                Used for helping the algorithm guessing the new position (default: {0})

        Raises:
            Exception: if the particle is outside the mesh

        Returns:
            [int] -- voronoi region index where the particle lies
        """
        position = self.findParticleVoronoiPosition(particle, lastKnownVoronoiPosition)

        if position is None:
            raise Exception(f"The particle is outside the mesh! {particle.x1}")

        # find region vertices coordinates
        coordinates = []
        for i in self.mesh.voronoiRegions[position]:
            coordinates.append(self.mesh.voronoi.vertices[i])

        # deposit particle to voronoi mesh
        w = computeBarycentric(particle.x1, coordinates)
        for i, weight in enumerate(w):
            # important: f lives in the primary mesh (i.e. triangles = simplices)
            # but for the python voronoi library ->
            # index(voronoi_vertex) = index(delaunay_simplex)
            vertexIndex = self.mesh.voronoiRegions[position][i]

            # find the area of the simplex, the weight must be normalized with the areas
            area = self.mesh.areaSimplex(vertexIndex)

            # update the distribution function
            self.f[vertexIndex] = self.f[vertexIndex] + weight / area

        return position

    def findParticleVoronoiPosition(self, particle: Particle,
                                    lastKnownVoronoiPosition=0):
        """Find the voronoi region in which the particle lies

        Arguments:
            particle {Particle} -- Particle object

        Keyword Arguments:
            lastKnownVoronoiPosition {int} -- last known voronoi position of
                the particle (default: {0})
        Returns:
            int -- Index of the voronoi regions in which the particle
        """
        checkedRegions: list(int) = []

        # start searching from the particle's last known position
        regionsToCheck: list(int) = [lastKnownVoronoiPosition]

        position: int = None
        while (len(checkedRegions) < len(self.mesh.voronoiRegions)
               and position is None):

            for i in regionsToCheck:
                if (i in checkedRegions):
                    continue
                elif (self.mesh.isPointInsideVoronoiRegion(particle.x1, i)):
                    position = i
                    break
                checkedRegions.append(i)

            # extend the search to the closest neighbors
            regionsToCheck = self.mesh.findVoronoiRegionsNeighbors(checkedRegions)

        return position
