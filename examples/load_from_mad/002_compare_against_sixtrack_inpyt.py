import pickle

import pysixtrack
import sixtracktools

with open('line_from_mad.pkl', 'rb') as fid:
    line_from_mad = pysixtrack.Line.from_dict(pickle.load(fid))

sixinput = sixtracktools.sixinput.SixInput('sixtrack/')
line_from_sixtrack, other = pysixtrack.Line.from_sixinput(sixinput)

bb_ele_mad, bb_names_mad = line_from_mad.get_elements_of_type(
        [pysixtrack.elements.BeamBeam4D, pysixtrack.elements.BeamBeam6D])
bb_ele_six, bb_names_six = line_from_sixtrack.get_elements_of_type(
        [pysixtrack.elements.BeamBeam4D, pysixtrack.elements.BeamBeam6D])


assert(len(bb_ele_six) == len(bb_ele_mad) == len(bb_names_mad) == len(bb_names_six))

for ee_mad, ee_six, nn_mad, nn_six in zip(bb_ele_mad, bb_ele_six, 
        bb_names_mad, bb_names_six):
    assert(nn_mad == nn_six)
    if nn_mad == 'bb_par.r1b1_5':
        prrrrr




