import numpy as np
from singleParticleIntegrator import Particle
# import matplotlib.pyplot as plt
# from scipy.spatial import voronoi_plot_2d
# from triangular_mesh import generateMesh

k = 1.0
m = 1.0
nsteps = 1
h = 1E+2
# nparticles = 10

# initial condition of the particle
# FIXME: generalize to n particles
x0 = np.array([0, 0])
x1 = np.array([2, 1])

# inizialize the particle
particle = Particle(x0, x1, m, k, h)

out = open("out.txt", "w+")

# compute initial energy
E0 = particle.energy()

for i in range(nsteps):
    # evolve the particle in time
    particle.particlePush()

    # energy error
    dE = (particle.energy() - E0) / E0

    out.write(str(i) + ' ' + str(dE) + " ")
    np.savetxt(out, particle.x1, newline=' ')
    out.write("\n")

# mesh = generateMesh(5, 10, 10)
# plt.triplot(mesh.delaunay.points[:, 0], mesh.delaunay.points[:, 1],
#             mesh.delaunay.simplices.copy())
# plt.plot(mesh.delaunay.points[:, 0], mesh.delaunay.points[:, 1], 'o')
# voronoi_plot_2d(mesh.voronoi)
# plt.show()
