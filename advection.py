#!/usr/bin/python3
import numpy as np
from singleParticleIntegrator import Particle
# import matplotlib.pyplot as plt
# from scipy.spatial import voronoi_plot_2d
from mesh import EquilateralTriangularMesh
from particleMeshDeposit import ParticlesMeshDeposit
import random

k = 1.0
m = 1.0
nsteps = 1000
h = 1E+1
nparticles = 10

particles = []
lastParticlesPosition = [0]*nparticles
for i in range(nparticles):
    # initial condition of the particle
    x0 = np.array([random.uniform(-2, 2), random.uniform(-2, 2)])
    x1 = np.array([random.uniform(-2, 2), random.uniform(-2, 2)])

    particle = Particle(x0, x1, m, k, h)
    particles.append(particle)

out = open("out.txt", "w+")

print("Inizializing mesh...")
mesh = EquilateralTriangularMesh(-10, 10, -10, 10, 1)
print("===============")
print("Mesh size: " + str(len(mesh.delaunay.simplices)) + " triangles")
print("N Particles " + str(nparticles))
print("N timesteps " + str(nsteps))
print("===============\n")

for i in range(nsteps):
    print("Timestep " + str(i))

    # rebuild the density function at each time step
    particleDeposit = ParticlesMeshDeposit(mesh)

    # evolve the particles in time
    for i in range(nparticles):
        particle = particles[i]
        particle.particlePush()

        lastParticlesPosition[i] = particleDeposit.addParticle(particle, lastParticlesPosition[i])

    # WRITE TO FILE (TRAJECTORY)
    # out.write(str(i) + ' ' + str(dE) + " ")
    # np.savetxt(out, particle.x1, newline=' ')
    # out.write("\n")

# PLOT MESH
# plt.triplot(mesh.delaunay.points[:, 0], mesh.delaunay.points[:, 1],
#             mesh.delaunay.simplices.copy())
# plt.plot(mesh.delaunay.points[:, 0], mesh.delaunay.points[:, 1], 'o')
# voronoi_plot_2d(mesh.voronoi)
# plt.show()
