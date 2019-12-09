import numpy as np

import pysixtrack
import pysixtrack.be_beamfields

element_list = [
    pysixtrack.elements.Drift,
    pysixtrack.elements.DriftExact,
    pysixtrack.elements.Multipole,
    pysixtrack.elements.Cavity,
    pysixtrack.elements.SawtoothCavity,
    pysixtrack.elements.XYShift,
    pysixtrack.elements.SRotation,
    pysixtrack.elements.RFMultipole,
    pysixtrack.elements.BeamMonitor,
    pysixtrack.elements.DipoleEdge,
    pysixtrack.elements.Line,
    pysixtrack.elements.LimitRect,
    pysixtrack.elements.LimitEllipse,
    pysixtrack.elements.SpaceChargeCoasting,
    pysixtrack.elements.SpaceChargeBunched,
]


def test_track_all():
    for el in element_list:
        p = pysixtrack.Particles()
        el().track(p)


def test_track_spacecharge():
    x_co = 0.1
    y_co = -0.5
    sigma_x = 0.5
    sigma_y = 0.1
    el1 = pysixtrack.elements.SpaceChargeBunched(
        number_of_particles=1e11,
        bunchlength_rms=0.22,
        sigma_x=sigma_x,
        sigma_y=sigma_y,
        length=2.0,
        x_co=x_co,
        y_co=y_co,
    )
    line_density = el1.number_of_particles / (
        el1.bunchlength_rms * np.sqrt(2 * np.pi)
    )
    el2 = pysixtrack.elements.SpaceChargeCoasting(
        line_density=line_density,
        sigma_x=el1.sigma_x,
        sigma_y=el1.sigma_y,
        length=el1.length,
        x_co=el1.x_co,
        y_co=el1.y_co,
    )

    # test absolute kick for sigma_x > sigma_y
    x_offset = 0.2
    y_offset = 0.5
    p1 = pysixtrack.Particles()
    p1.x = x_co + x_offset
    p1.y = y_co + y_offset
    p2 = p1.copy()
    el1.track(p1)
    el2.track(p2)
    # print(p1.px,p1.py)
    assert np.isclose(p1.px, 1.33671170365656e-07, atol=1e-15) == True
    assert np.isclose(p1.py, 6.228154616877465e-07, atol=1e-15) == True
    assert p1.compare(p2, abs_tol=1e-15)

    el1.sigma_x = sigma_y
    el1.sigma_y = sigma_x
    el2.sigma_x = sigma_y
    el2.sigma_y = sigma_x
    p1 = pysixtrack.Particles()
    p1.x = x_co + y_offset
    p1.y = y_co + x_offset
    p2 = p1.copy()
    el1.track(p1)
    el2.track(p2)
    assert np.isclose(p1.px, 6.228154616877465e-07, atol=1e-15) == True
    assert np.isclose(p1.py, 1.33671170365656e-07, atol=1e-15) == True
    assert p1.compare(p2, abs_tol=1e-15)

    p1 = pysixtrack.Particles()
    p1.x = el1.x_co
    p1.y = el1.y_co
    p2 = p1.copy()
    el1.track(p1)
    el2.track(p2)
    assert np.isclose(p1.px, 0.0, atol=1e-15) == True
    assert np.isclose(p1.py, 0.0, atol=1e-15) == True
    assert p1.compare(p2, abs_tol=1e-15)


def test_track_rfmultipole():
    p1 = pysixtrack.Particles()
    p1.x = 1
    p1.y = 1
    p2 = p1.copy()

    el1 = pysixtrack.elements.RFMultipole(knl=[0.5, 2, 0.2], ksl=[0.5, 3, 0.1])
    el2 = pysixtrack.elements.Multipole(knl=el1.knl, ksl=el1.ksl)

    el1.track(p1)
    el2.track(p2)

    assert p1.compare(p2, abs_tol=1e-15)


def test_track_LimitRect():
    min_x = -0.1
    max_x = 0.3
    min_y = -0.5
    max_y = 0.1
    el = pysixtrack.elements.LimitRect(
        min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y
    )

    p1 = pysixtrack.Particles()
    p1.x = 1
    p1.y = 1
    ret = el.track(p1)
    assert ret == "Particle lost"

    arr = np.arange(0, 1, 0.001)
    p2 = pysixtrack.Particles(x=arr, y=arr)
    ret = el.track(p2)

    p2.x += max_x + 1e-6
    ret = el.track(p2)
    assert ret == "All particles lost"


def test_track_LimitEllipse():
    limit_a = 0.1
    limit_b = 0.2
    el = pysixtrack.elements.LimitEllipse(a=limit_a, b=limit_b)

    p1 = pysixtrack.Particles()
    p1.x = 1
    p1.y = 1
    ret = el.track(p1)
    assert ret == "Particle lost"

    arr = np.arange(0, 1, 0.001)
    p2 = pysixtrack.Particles(x=arr, y=arr)
    ret = el.track(p2)
    survived = np.where(
        (p2.x ** 2 / limit_a ** 2 + p2.y ** 2 / limit_b ** 2 <= 1.0)
    )
    assert len(p2.state) == len(survived[0])

    p2.x += limit_a + 1e-6
    ret = el.track(p2)
    assert ret == "All particles lost"
