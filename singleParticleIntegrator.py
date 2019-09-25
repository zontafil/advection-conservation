import numpy as np
class Particle:
    # initialize a particle with initial condition x0, x1,
    # mass m, timestep h
    def __init__(self, x0, x1, m, k, h):
        self.x0 = x0
        self.x1 = x1
        self.m = m
        self.h = h
        self.k = k

    def particlePush(self):
        # evolve particle using a variational midpoint rule symplectic integrator
        ap = self.m / self.h + self.h * self.k /4.0
        am = self.m / self.h - self.h * self.k /4.0

        xt = self.x1
        self.x1 =  (2.0 * self.x1 * am / ap) - self.x0
        self.x0 = xt

    def discreteLeftLegendre(self):
        return self.m / self.h * (self.x1 - self.x0) + self.h * self.k / 4.0 * (self.x0 + self.x1)

    def energy(self):
        p0 = self.discreteLeftLegendre()
        return (np.linalg.norm(p0)**2 /2.0/ self.m + 0.5 * self.k * np.linalg.norm(self.x0)**2)