import pickle
import numpy as np

import pysixtrack
import sixtracktools

with open('line_from_mad.pkl', 'rb') as fid:
    lmad = pysixtrack.Line.from_dict(pickle.load(fid))

sixinput = sixtracktools.sixinput.SixInput('sixtrack/')
lsix, other = pysixtrack.Line.from_sixinput(sixinput)

assert((lmad.get_length() - lsix.get_length())<1e-6)

for ll in (lmad, lsix):
    ll.remove_inactive_multipoles(inplace=True)
    ll.remove_zero_length_drifts(inplace=True)
    ll.merge_consecutive_drifts(inplace=True)

for ii, (eem, ees) in enumerate(zip(lmad.elements, lsix.elements)):
    assert(type(eem) == type(ees))
