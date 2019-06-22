import numpy as np
import pickle

from cpymad.madx import Madx
import pysixtrack

from pysixtrack.be_beambeam.tools import norm, find_alpha_and_phi
from pysixtrack.be_beambeam.tools import get_bb_names_madpoints_sigmas
from pysixtrack.be_beambeam.tools import shift_strong_beam_based_on_close_ip 
from pysixtrack.be_beambeam.tools import setup_beam_beam_in_line 
from pysixtrack import MadPoint

mad=Madx()
mad.options.echo=False;
mad.options.warn=False;
mad.options.info=False;


# Load sequence
mad.call('mad/lhcwbb.seq')

# Parameters to be cross-checked
n_slices = 11
bb_coupling = False

# Disable beam-beam kicks 
mad.globals.on_bb_charge = 0.

ip_names = [1, 2, 5, 8]

# Retrieve geometry information

# IP locations from the survey
mad.use('lhcb1'); mad.twiss(); mad.survey()
IP_xyz_b1 = {}
for ip in ip_names:
    IP_xyz_b1[ip] = MadPoint('ip%d'%ip+':1', mad, add_CO=False)

mad.use('lhcb2'); mad.twiss(); mad.survey()
IP_xyz_b2 = {}
for ip in ip_names:
    IP_xyz_b2[ip] = MadPoint('ip%d'%ip+':1', mad, add_CO=False)

# Beam-beam names and locations
bb_names_b1, bb_xyz_b1, bb_sigmas_b1 = get_bb_names_madpoints_sigmas(
        mad, seq_name='lhcb1')
bb_names_b2, bb_xyz_b2, bb_sigmas_b2 = get_bb_names_madpoints_sigmas(
        mad, seq_name='lhcb2')

# Get bunch paramenters
bunch_intensity = mad.sequence.lhcb1.beam.npart
gamma_r = mad.sequence.lhcb1.beam.gamma
beta_r = np.sqrt(1-1./gamma_r**2)


# Check naming convention
assert len(bb_names_b1)==len(bb_names_b2)
for nbb1, nbb2 in zip(bb_names_b1, bb_names_b2):
    assert(nbb1==nbb2.replace('b2_','b1_'))

# Check number of slices
assert(len(
    [nn for nn in bb_names_b1 if nn.startswith('bb_ho.l1')])==(n_slices-1)/2)

shift_strong_beam_based_on_close_ip(
        points_weak=bb_xyz_b1, 
        points_strong=bb_xyz_b2,
        IPs_survey_weak=IP_xyz_b1,
        IPs_survey_strong=IP_xyz_b2)

line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.lhcb1)

mad_ft = Madx()
mad_ft.options.echo=False;
mad_ft.options.warn=False;
mad_ft.options.info=False;

mad_ft.call('mad/lhcwbb_fortracking.seq')
# without this the sequence does not work properly
mad_ft.use('lhcb1')

line_for_tracking, _ = pysixtrack.Line.from_madx_sequence(
        mad_ft.sequence['lhcb1'])

setup_beam_beam_in_line(line, bb_names_b1,
        bb_sigmas_strong=bb_sigmas_b2,
        bb_points_weak=bb_xyz_b1, bb_points_strong=bb_xyz_b2,
        bunch_intensity=bunch_intensity, n_slices_6D=n_slices)

ready_bb_elems = line.get_elements_of_type(
        (pysixtrack.elements.BeamBeam4D, 
         pysixtrack.elements.BeamBeam6D))
dct_ready_bb = {nn:ee for nn, ee in zip(ready_bb_elems[1], ready_bb_elems[0])}

for ii, nn in enumerate(line_for_tracking.element_names):
    if nn in dct_ready_bb.keys():
        line_for_tracking.elements[ii] = dct_ready_bb[nn]

assert(np.abs(line_for_tracking.get_length()\
        - mad.sequence.lhcb1.beam.circ)<1e-6)

# There is a problem in the mask 
# (the RF frequancy is wrong in the machine for tracking
# I patch it here -> to be fixed properly!!!!
dct_correct_cavities = dict(zip(*line.get_elements_of_type(
    pysixtrack.elements.Cavity)[::-1]))
for ii, nn in enumerate(line_for_tracking.element_names):
    if nn in dct_correct_cavities.keys():
        line_for_tracking.elements[ii].frequency = \
                dct_correct_cavities[nn].frequency

with open('line_from_mad.pkl', 'wb') as fid:
    pickle.dump(line_for_tracking.to_dict(keepextra=True), fid)

import matplotlib.pyplot as plt
plt.close('all')
fig1 = plt.figure(1)

plt.plot(           [pb.p[2] for pb in bb_xyz_b1],
                    [pb.p[0] for pb in bb_xyz_b1], 'b.')
plt.quiver(np.array([pb.p[2] for pb in bb_xyz_b1]),
           np.array([pb.p[0] for pb in bb_xyz_b1]), 
          np.array([pb.ex[2] for pb in bb_xyz_b1]),
          np.array([pb.ex[0] for pb in bb_xyz_b1]))
plt.quiver(np.array([pb.p[2] for pb in bb_xyz_b1]),
           np.array([pb.p[0] for pb in bb_xyz_b1]), 
          np.array([pb.ez[2] for pb in bb_xyz_b1]),
          np.array([pb.ez[0] for pb in bb_xyz_b1]))

plt.plot(           [pb.p[2] for pb in bb_xyz_b2],
                    [pb.p[0] for pb in bb_xyz_b2], 'r.')
plt.quiver(np.array([pb.p[2] for pb in bb_xyz_b2]),
           np.array([pb.p[0] for pb in bb_xyz_b2]), 
          np.array([pb.ex[2] for pb in bb_xyz_b2]),
          np.array([pb.ex[0] for pb in bb_xyz_b2]))
plt.quiver(np.array([pb.p[2] for pb in bb_xyz_b2]),
           np.array([pb.p[0] for pb in bb_xyz_b2]), 
          np.array([pb.ez[2] for pb in bb_xyz_b2]),
          np.array([pb.ez[0] for pb in bb_xyz_b2]))
plt.show()


