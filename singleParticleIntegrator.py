import numpy as np


class Particle:
    """Particle object with 2-dim harmonic oscillator
    symplectic integrator included

    Attributes:
        x0 {numpy.array} -- position at time step k
        x1 {numpy.array} -- position at time step k + 1
        m {int} -- mass of the particle
        k {int} -- elastic constant
        h {int} -- time step
    """
    def __init__(self, x0, x1, m, k, h):
        """Initialize the particle

        Arguments:
            x0 {numpy.array} -- position at time step k
            x1 {numpy.array} -- position at time step k + 1
            m {int} -- mass of the particle
            k {int} -- elastic constant
            h {int} -- time step
        """
        self.x0 = x0
        self.x1 = x1
        self.m = m
        self.h = h
        self.k = k

    def particlePush(self):
        """Evolve particle using a variational midpoint rule symplectic integrator
        """
        ap = self.m / self.h + self.h * self.k / 4.0
        am = self.m / self.h - self.h * self.k / 4.0

        xt = self.x1
        self.x1 = (2.0 * self.x1 * am / ap) - self.x0
        self.x0 = xt

    def discreteLeftLegendre(self):
        """Compute the discrete left Legendre transform

        Returns:
            numpy.array -- discrete conjugate momenta at time step k (p0)
        """
        return self.m / self.h * (self.x1 - self.x0) \
            + self.h * self.k / 4.0 * (self.x0 + self.x1)

    def energy(self):
        """compute the energy of the particle at x0,p0

        Returns:
            int -- energy of the particle
        """
        p0 = self.discreteLeftLegendre()
        return (np.linalg.norm(p0)**2 / 2.0 / self.m
                + 0.5 * self.k * np.linalg.norm(self.x0)**2)
