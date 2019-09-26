from mesh import Mesh
from singleParticleIntegrator import Particle


class ParticlesMeshDeposit:
    """Deposit and distribution function matrices and methods

    Attributes:
        mesh
        f -- distribution function
        voronoiWeights -- Weights of the particles in the dual mesh

    """
    def __init__(self, mesh: Mesh):
        """Constructor. A valid mesh must be passed

        Arguments:
            mesh {Mesh}
        """
        self.mesh = mesh
        self.f: list(float) = [None]*len(mesh.delaunay.points)
        self.voronoiWeights: list(float) = [None]*len(mesh.voronoi.vertices)

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

        # TODO

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
