import pickle
import numpy as np

import pysixtrack
import sixtracktools

with open('line_from_mad.pkl', 'rb') as fid:
    lmad = pysixtrack.Line.from_dict(pickle.load(fid))

sixinput = sixtracktools.sixinput.SixInput('sixtrack/')
lsix, other = pysixtrack.Line.from_sixinput(sixinput)

bucket_length = lmad.get_length()/35640
lr_spacing = 5*bucket_length

assert((lmad.get_length() - lsix.get_length())<1e-6)

i_bb = 1
s_ip5_six = lsix.get_s_elements()[lsix.element_names.index('ip5')]
s_bb_six = lsix.get_s_elements()[lsix.element_names.index('bb_par.l5b1_%d'%i_bb)]
nsix = (s_bb_six-s_ip5_six)/lr_spacing 
s_ip5_mad = lmad.get_s_elements()[lmad.element_names.index('ip5')]
s_bb_mad = lmad.get_s_elements()[lmad.element_names.index('bb_par.l5b1_%d'%i_bb)]
nmad = (s_bb_mad-s_ip5_mad)/lr_spacing 

print(f"nmad={nmad} six={nsix}")

# for ll in (lmad, lsix):
#     ll.remove_inactive_multipoles(inplace=True)
#     ll.remove_zero_length_drifts(inplace=True)
#     ll.merge_consecutive_drifts(inplace=True)
# 
# for ii, (eem, ees) in enumerate(zip(lmad.elements, lsix.elements)):
#     assert(type(eem) == type(ees))
