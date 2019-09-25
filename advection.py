import numpy as np
from singleParticleIntegrator import Particle
# import scipy
# import matplotlib.pyplot as plt
# from triangular_mesh import generateMesh
# from scipy.spatial import Delaunay
# from scipy.spatial import Voronoi, voronoi_plot_2d

k = 1.0
m = 1.0
nsteps = 1
h = 1E+2
# nparticles = 10

# initial condition of the particle
# FIXME: generalize to n particles
x0 = np.array([0,0])
x1 = np.array([2,1])

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
    
# points = generateMesh(5, 10, 10)
# tri = Delaunay(points)
# plt.triplot(points[:,0], points[:,1], tri.simplices.copy())
# plt.plot(points[:,0], points[:,1], 'o')
# vor = Voronoi(points)
# voronoi_plot_2d(vor)
# plt.show()
