import numpy as np

from cpymad.madx import Madx
import pysixtrack
from pysixtrack.particles import Particles

# run MADX tests
mad = Madx()
mad.call("rf_multipole.madx")

# same test in pysixtrack
mad_sequence = mad.sequence["sequ_rfmultipole"]
rf_mulitpole_mad = mad_sequence.elements[1]
freq = rf_mulitpole_mad.freq * 1e6  # MAD units are MHz
knl = rf_mulitpole_mad.knl
pn = np.array(rf_mulitpole_mad.pnl) * 360  # MAD units are 2pi
lag = rf_mulitpole_mad.lag * 360  # MAD units are 2pi

rf_multipole = pysixtrack.elements.RFMultipole(
    voltage=0, frequency=freq, lag=lag, knl=knl, ksl=[0], pn=pn, ps=[0]
)

mad_part=Particles.from_mad_track(mad)
p1=mad_part.copy(0)
p2=mad_part.copy(1)
p3=p1.copy()
rf_multipole.track(p3)

p2.compare(p3)


