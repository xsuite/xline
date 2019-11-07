import pysixtrack

element_list = [
    pysixtrack.elements.Drift,
    pysixtrack.elements.DriftExact,
    pysixtrack.elements.Multipole,
    pysixtrack.elements.Cavity,
    pysixtrack.elements.XYShift,
    pysixtrack.elements.SRotation,
    pysixtrack.elements.RFMultipole,
    pysixtrack.elements.BeamMonitor,
    pysixtrack.elements.DipoleEdge,
    pysixtrack.elements.Line,
    pysixtrack.elements.LimitRect,
    pysixtrack.elements.LimitEllipse
]



def test_track_all():
    for el in element_list:
        p = pysixtrack.Particles()
        el().track(p)

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
