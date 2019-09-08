import numpy as np
from pysixtrack.particles import Particles

p = Particles()
p.x = np.arange(10, dtype=np.float64)

# I cut away all the odd coordinates (to check that I can skip)
p.state = np.int_(np.mod(np.int_(p.x), 2) == 0)
p.remove_lost_particles()

# I cut away internal particles
p.state = np.int_(p.x>5)
p.remove_lost_particles()

assert(np.all(p.x == np.array([6., 8.])))
assert(np.all(p.lost_particles[0].x == np.array([1,3,5,7,9])))
assert(np.all(p.lost_particles[1].x == np.array([0, 2.,4.])))

