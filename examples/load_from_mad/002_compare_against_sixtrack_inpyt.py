import pickle
import numpy as np

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

    dmad = ee_mad.to_dict(keepextra=True)
    dsix = ee_six.to_dict(keepextra=True)

    for kk in dmad.keys():
        if dmad[kk] == dsix[kk]:
            continue

        try:
            diff_rel = np.abs((dmad[kk] - dsix[kk])/dmad[kk])
        except ZeroDivisionError:
            diff_rel = 100.
        check_rel = diff_rel < 1e-5
        if not check_rel:
            diff_abs = np.abs(dmad[kk] - dsix[kk])
            check_abs = diff_abs< 1e-12
            print(f"{nn_mad}[{kk}] - mad:{dmad[kk]} six:{dsix[kk]}")

            # Exceptions:
            if kk == 'x_bb':
                if diff_abs/dmad['sigma_x'] < 0.001:
                    continue
            if kk == 'y_bb':
                if diff_abs/dmad['sigma_y'] < 0.001:
                    continue
            if kk=='alpha':
                if diff_abs<10e-6:
                    continue
            if kk=='x_bb_co':
                if diff_abs/np.sqrt(dmad['sigma_11']) < 0.001:
                    continue
            if kk=='y_bb_co':
                if diff_abs/np.sqrt(dmad['sigma_33']) < 0.001:
                    continue

            # general check
            if not check_abs:
                raise ValueError('Too large discrepancy!')


