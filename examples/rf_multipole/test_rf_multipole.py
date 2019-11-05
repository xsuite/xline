import numpy as np

from cpymad.madx import Madx
import pysixtrack
from pysixtrack.particles import Particles

# run MADX tests
mad = Madx()
mad.call("test_rf_multipole.madx")
tracksumm_madx = mad.table.tracksumm

# same test in pysixtrack
mad_sequence = mad.sequence["sequ_rfmultipole"]
rf_mulitpole_mad = mad_sequence.elements[1]
freq = rf_mulitpole_mad.freq * 1e6  # MAD units are MHz
knl = rf_mulitpole_mad.knl
pnl = np.array(rf_mulitpole_mad.pnl) * 360  # MAD units are 2pi
lag = rf_mulitpole_mad.lag * 360  # MAD units are 2pi

my_rf_multipole = pysixtrack.elements.RFMultipole(
    voltage=0, frequency=freq, lag=lag, knl=knl, ksl=[0], pn=pnl, ps=[0]
)

p0c = mad_sequence.beam.pc * 1e9
x = tracksumm_madx.x[0]
px = tracksumm_madx.px[0]
y = tracksumm_madx.y[0]
py = tracksumm_madx.py[0]
t = tracksumm_madx.t[0]
pt = tracksumm_madx.pt[0]

part = Particles(p0c=p0c, x=x, px=px, y=y, py=py, tau=t, pt=pt)
# print(part)

my_rf_multipole.track(part)
print(part)

# part(p0c=p0c, x=x, px=px, y=y, py=py, zeta=?, delta=?)
""
