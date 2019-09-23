import numpy as np

k = 1.0
m = 1.0
nsteps = 10000
h = 1E+2

# initial condition of the particle
x0 = np.array([0,0])
x1 = np.array([2,1])

def particlePush(x0, x1, h):
    # evolve particle using a variational midpoint rule symplectic integrator
    ap = m / h + h*k/4.0
    am = m / h - h*k/4.0

    return (2.0 * x1 * am / ap) - x0

def discreteLeftLegendre(x0, x1, h):
    return m / h * (x1 - x0) + h * k / 4 * (x0 + x1)

def energy(x0, x1, h):
    p0 = discreteLeftLegendre(x0, x1, h)
    return (np.linalg.norm(p0)**2 /2.0/m + 0.5 * k * np.linalg.norm(x0)**2)

out = open("out.txt", "w+")

# compute initial energy
E0 = energy(x0, x1, h)

for i in range(nsteps):
    xt = x1
    x1 = particlePush(x0, x1, h)
    x0 = xt

    E = energy(x0, x1, h)
    # energy error
    dE = (energy(x0, x1, h) - E0) / E0

    # out.write(np.array_str(x1))
    out.write(str(i) + ' ' + str(dE) + " ")
    np.savetxt(out, x1, newline=' ')
    out.write("\n")
    

