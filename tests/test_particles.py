from xline import XlineTestParticles


def test_particle_init():
    p = XlineTestParticles()
    print(p)


def test_particle_reference():
    p = XlineTestParticles(p0c=1e9)
    p.beta0 = 0.91
    err = abs(p.p0c ** 2 + p.mass0 ** 2 - p.energy0 ** 2) / p.mass0 ** 2
    assert err < 1e-15
    p.beta0 = 0.9101
    err = abs(p.p0c ** 2 + p.mass0 ** 2 - p.energy0 ** 2) / p.mass0 ** 2
    assert err < 1.17e-15
    p.gamma0 = 1.99
    err = abs(p.p0c ** 2 + p.mass0 ** 2 - p.energy0 ** 2) / p.mass0 ** 2
    assert err < 1e-15
    # p.energy0 = 2 * p.mass0
    # err = abs(p.p0c ** 2 + p.mass0 ** 2 - p.energy0 ** 2) / p.mass0 ** 2
    # assert err == 0
    p.p0c = 0.1 * p.mass0
    err = abs(p.p0c ** 2 + p.mass0 ** 2 - p.energy0 ** 2) / p.mass0 ** 2
    assert err < 1e-15
