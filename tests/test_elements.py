from pysixtrack import elements


def test_element_creation():
   el = elements.Drift()
   el = elements.DriftExact()
   el = elements.Multipole()
   assert el.order == 0
   el = elements.XYShift()
   el = elements.SRotation()
   el = elements.Cavity()
   el = elements.Line()
   el = elements.DipoleEdge()
