import numpy as np
from xline.particles import XlineTestParticles


def test_particle_loss():
    p = XlineTestParticles(p0c=1e9,
            x=np.arange(10, dtype=np.float64))

    # I cut away all the odd coordinates (to check that I can skip)
    p.state = np.int_(np.mod(np.int_(p.x), 2) == 0)
    p.remove_lost_particles()

    # I cut away internal particles
    p.state = np.int_(p.x > 5)
    p.remove_lost_particles()

    assert np.all(p.x == np.array([6.0, 8.0]))
