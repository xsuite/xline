import numpy as np

from cpymad.madx import Madx
import pysixtrack

# run MADX tests
mad = Madx()
mad.call("rf_multipole.madx")

# create pysixtrack rfmultipole
mad_sequence = mad.sequence["sequ_rfmultipole"]
rf_mulitpole_mad = mad_sequence.elements[1]
freq = rf_mulitpole_mad.freq * 1e6  # MAD units are MHz
knl = rf_mulitpole_mad.knl
pn = np.array(rf_mulitpole_mad.pnl) * 360  # MAD units are 2pi
lag = rf_mulitpole_mad.lag * 360  # MAD units are 2pi

rf_multipole = pysixtrack.elements.RFMultipole(
    voltage=0, frequency=freq, lag=lag, knl=knl, ksl=[0], pn=pn, ps=[0]
)


# track pysixtrack
mad_part = pysixtrack.Particles.from_madx_track(mad)
p1 = mad_part.copy(0)
p2 = mad_part.copy(1)
p3 = p1.copy()
rf_multipole.track(p3)

# compare
p2.compare(p3)

# test conversion
line = pysixtrack.Line.from_madx_sequence(mad_sequence)
tw = mad.twiss(betx=1, bety=1, x=0.1, t=0.5)
p_mad = pysixtrack.Particles.from_madx_twiss(tw)
p_six = mad_part.copy(0)
p_out = pysixtrack.Particles.from_list(
    line.track_elem_by_elem(p_six, start=False)
)
