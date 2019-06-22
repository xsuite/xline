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

assert(len(lmad) == len(lsix))

assert((lmad.get_length() - lsix.get_length())<1e-6)

def norm(x):
    return np.sqrt(np.sum(np.array(x)**2))

for ii, (ee_mad, ee_six, nn_mad, nn_six) in enumerate(zip(
    lmad.elements, lsix.elements, lmad.element_names, lsix.element_names)):
    assert(type(ee_mad) == type(ee_six))

    dmad = ee_mad.to_dict(keepextra=True)
    dsix = ee_six.to_dict(keepextra=True)

    for kk in dmad.keys():
        
        # Check if they are identical
        if dmad[kk] == dsix[kk]:
            continue
        
        # Check if the relative error is small
        try:
            diff_rel = norm(
                    np.array(dmad[kk]) - np.array(dsix[kk]))/norm(dmad[kk])
        except ZeroDivisionError:
            diff_rel = 100.
        if diff_rel < 3e-5: 
            continue

        # Check if absolute error is small
        diff_abs = norm(np.array(dmad[kk]) - np.array(dsix[kk]))
        print(f"{nn_mad}[{kk}] - mad:{dmad[kk]} six:{dsix[kk]}")
        if diff_abs < 1e-12:
            continue

        # Exception: multipole lrad is not passes to sixtraxk
        if isinstance(ee_mad, pysixtrack.elements.Multipole):
            if kk == 'length':
                continue


            ## Exceptions:
            #if kk == 'x_bb':
            #    if diff_abs/dmad['sigma_x'] < 0.001:
            #        continue
            #if kk == 'y_bb':
            #    if diff_abs/dmad['sigma_y'] < 0.001:
            #        continue
            #if kk=='alpha':
            #    if diff_abs<10e-6:
            #        continue
            #if kk=='x_bb_co':
            #    if diff_abs/np.sqrt(dmad['sigma_11']) < 0.001:
            #        continue
            #if kk=='y_bb_co':
            #    if diff_abs/np.sqrt(dmad['sigma_33']) < 0.001:
            #        continue

            # general check
    
        # If it got here it means that no condition above is met
        raise ValueError('Too large discrepancy!')

