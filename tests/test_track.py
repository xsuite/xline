import numpy as np

import xline

element_list = [
    xline.elements.Drift,
    xline.elements.DriftExact,
    xline.elements.Multipole,
    xline.elements.Cavity,
    xline.elements.SawtoothCavity,
    xline.elements.XYShift,
    xline.elements.SRotation,
    xline.elements.RFMultipole,
    xline.elements.BeamMonitor,
    xline.elements.DipoleEdge,
    xline.elements.Line,
    xline.elements.LimitRect,
    xline.elements.LimitEllipse,
    xline.elements.LimitRectEllipse,
    xline.elements.BeamBeam4D,
    xline.elements.BeamBeam6D,
    xline.elements.SCCoasting,
    xline.elements.SCQGaussProfile,
]


def test_track_all():
    for el in element_list:
        p = xline.XlineTestParticles(p0c=1e9)
        el().track(p)


def test_track_rfmultipole():
    p1 = xline.XlineTestParticles(p0c=1e9)
    p1.x = 1
    p1.y = 1
    p2 = p1.copy()

    el1 = xline.elements.RFMultipole(knl=[0.5, 2, 0.2], ksl=[0.5, 3, 0.1])
    el2 = xline.elements.Multipole(knl=el1.knl, ksl=el1.ksl)

    el1.track(p1)
    el2.track(p2)

    assert p1.compare(p2, abs_tol=1e-15)


def test_track_LimitRect():
    min_x = -0.1
    max_x = 0.3
    min_y = -0.5
    max_y = 0.1
    el = xline.elements.LimitRect(
        min_x=min_x, max_x=max_x, min_y=min_y, max_y=max_y
    )

    p1 = xline.XlineTestParticles(p0c=1e9)
    p1.x = 1
    p1.y = 1
    ret = el.track(p1)
    assert p1.state == 0

    arr = np.arange(0, 1, 0.001)
    p2 = xline.XlineTestParticles(x=arr, y=arr)
    survive = np.where(
        (p2.x >= min_x) & (p2.x <= max_x) & (p2.y >= min_y) & (p2.y <= max_y)
    )
    ret = el.track(p2)
    assert len(p2.state) == len(survive[0])

    p2.x += max_x + 1e-6
    el.track(p2)
    assert len(p2.x) == 0


def test_track_LimitEllipse():
    limit_a = 0.1
    limit_b = 0.2
    el = xline.elements.LimitEllipse(a=limit_a, b=limit_b)

    p1 = xline.XlineTestParticles(p0c=1e9)
    p1.x = 1
    p1.y = 1
    ret = el.track(p1)
    assert p1.state == 0

    arr = np.arange(0, 1, 0.001)
    p2 = xline.XlineTestParticles(x=arr, y=arr)
    survive = np.where(
        (p2.x ** 2 / limit_a ** 2 + p2.y ** 2 / limit_b ** 2 <= 1.0)
    )
    ret = el.track(p2)
    assert len(p2.state) == len(survive[0])

    p2.x += limit_a + 1e-6
    ret = el.track(p2)
    assert len(p2.x) == 0


def test_track_LimitRectEllipse():
    limit_a = 0.1
    limit_b = 0.2
    max_x = 0.1
    max_y = 0.05
    el = xline.elements.LimitRectEllipse(
        max_x=max_x, max_y=max_y, a=limit_a, b=limit_b
    )

    p1 = xline.XlineTestParticles()
    p1.x = 1
    p1.y = 1
    ret = el.track(p1)
    assert p1.state == 0

    arr = np.arange(0, 1, 0.001)
    p2 = xline.XlineTestParticles(x=arr, y=arr)
    survive = np.where(
        (p2.x ** 2 / limit_a ** 2 + p2.y ** 2 / limit_b ** 2 <= 1.0)
        & (p2.x >= -max_x)
        & (p2.x <= max_x)
        & (p2.y >= -max_y)
        & (p2.y <= max_y)
    )
    ret = el.track(p2)
    assert len(p2.state) == len(survive[0])

    p2.x += limit_a + 1e-6
    ret = el.track(p2)
    assert len(p2.x) == 0
