import pickle
import numpy as np

import pysixtrack
import sixtracktools

with open('line_from_mad.pkl', 'rb') as fid:
    line_from_mad = pysixtrack.Line.from_dict(pickle.load(fid))

sixinput = sixtracktools.sixinput.SixInput('sixtrack/')
line_from_sixtrack, other = pysixtrack.Line.from_sixinput(sixinput)

assert((line_from_mad.get_length() - line_from_sixtrack.get_length())<1e-6)

for ll in (line_from_mad, line_from_sixtrack):
    ll.remove_inactive_multipoles(inplace=True)
    ll.remove_zero_length_drifts(inplace=True)
    ll.merge_consecutive_drifts(inplace=True)
