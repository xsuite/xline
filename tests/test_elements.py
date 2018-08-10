from pysixtrack import *


el = Drift()
el = DriftExact()
el = Multipole()
assert el.order == 0
el = XYShift()
el = SRotation()
el = Cavity()
el = Line()
