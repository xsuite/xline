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

i_bb = 21
s_ip5_six = lsix.get_s_elements()[lsix.element_names.index('ip5')]
s_bb_six = lsix.get_s_elements()[lsix.element_names.index('bb_par.l5b1_%d'%i_bb)]
nsix = (s_bb_six-s_ip5_six)/lr_spacing 
s_ip5_mad = lmad.get_s_elements()[lmad.element_names.index('ip5')]
s_bb_mad = lmad.get_s_elements()[lmad.element_names.index('bb_par.l5b1_%d'%i_bb)]
nmad = (s_bb_mad-s_ip5_mad)/lr_spacing 

print(f"nmad={nmad} six={nsix}")

import cpymad.madx
mad = cpymad.madx.Madx()
mad.call('mad/lhcwbb_fortracking.seq')
mad.use('lhcb1')# Without this the sequence is corrupted
testline, _ = pysixtrack.Line.from_madx_sequence(mad.sequence.lhcb1)
s_bb_test = testline.get_s_elements()[testline.element_names.index(
    'bb_par.l5b1_%d'%i_bb)]
print(s_bb_six, s_bb_test)

mad.use('lhcb1'); mad.twiss()
mad_names = list(mad.table.twiss.name)
mad_s = list(mad.table.twiss.s)

i_bb_tw = mad_names.index('bb_par.l5b1_%d:1'%i_bb)
s_bb_tw = mad_s[i_bb_tw]

seq = mad.sequence.lhcb1
s_seq = list(seq.element_positions())
n_seq = list(seq.element_names())
s_bb_seq = s_seq[n_seq.index('bb_par.l5b1_%d'%i_bb)] 

# for ll in (lmad, lsix):
#     ll.remove_inactive_multipoles(inplace=True)
#     ll.remove_zero_length_drifts(inplace=True)
#     ll.merge_consecutive_drifts(inplace=True)
# 
# for ii, (eem, ees) in enumerate(zip(lmad.elements, lsix.elements)):
#     assert(type(eem) == type(ees))
