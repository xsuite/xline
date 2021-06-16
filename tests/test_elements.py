from xline import elements


def test_other():
    el = elements.XYShift()
    el = elements.SRotation()
    el = elements.Cavity()
    el = elements.Line()
    el = elements.DipoleEdge()
    print(el)


def test_drift():
    el = elements.Drift(length=4)
    assert el.length == 4
    el = elements.DriftExact(length=4)
    assert el.length == 4


def test_multipole():
    el = elements.Multipole()
    assert el.order == 0
    assert el.knl[el.order] == 0
    assert el.ksl[el.order] == 0
    el = elements.Multipole(knl=[1])
    assert el.knl == [1]
    assert el.ksl == [0]
    assert el.order == 0
    el = elements.Multipole(knl=[1, 2, 3])
    assert el.order == 2
    el = elements.Multipole(ksl=[1, 2, 3])
    assert el.order == 2
