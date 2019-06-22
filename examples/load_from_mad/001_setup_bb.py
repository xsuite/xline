import numpy as np
import pickle

from cpymad.madx import Madx

import pysixtrack

from pysixtrack.be_beambeam.tools import norm, find_alpha_and_phi
from pysixtrack.be_beambeam.tools import get_bb_names_xyz_points_sigma_matrices
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
bb_names_b1, bb_xyz_b1, bb_sigmas_b1 = get_bb_names_xyz_points_sigma_matrices(
        mad, seq_name='lhcb1')
bb_names_b2, bb_xyz_b2, bb_sigmas_b2 = get_bb_names_xyz_points_sigma_matrices(
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

# Check that reference systems are parallel at the IPs
for ip in ip_names:
    assert(1. - np.dot(IP_xyz_b1[ip].ex, IP_xyz_b2[ip].ex) <1e-12)
    assert(1. - np.dot(IP_xyz_b1[ip].ey, IP_xyz_b2[ip].ey) <1e-12)
    assert(1. - np.dot(IP_xyz_b1[ip].ez, IP_xyz_b2[ip].ez) <1e-12)

# Shift B2 so that survey is head-on at the closest IP
# and find vector separation
sep_x = []
sep_y = []
for i_bb, name_bb in enumerate(bb_names_b1):
    
    pb1 = bb_xyz_b1[i_bb]
    pb2 = bb_xyz_b2[i_bb]
    
    # Find closest IP
    d_ip = 1e6
    use_ip = 0
    for ip in ip_names:
        dd = norm(pb1.p - IP_xyz_b1[ip].p)
        if dd < d_ip:
            use_ip = ip
            d_ip = dd

    # Shift B2
    shift_12 = IP_xyz_b2[use_ip].p - IP_xyz_b1[use_ip].p
    pb2.p -= shift_12

    # Find v12
    vbb_12 = bb_xyz_b2[i_bb].p - bb_xyz_b1[i_bb].p

    # Check that the two reference system are parallel
    try:
        assert(norm(pb1.ex-pb2.ex)<1e-10) #1e-4 is a reasonable limit
        assert(norm(pb1.ey-pb2.ey)<1e-10) #1e-4 is a reasonable limit
        assert(norm(pb1.ez-pb2.ez)<1e-10) #1e-4 is a reasonable limit
    except AssertionError:
        print(name_bb, 'Reference systems are not parallel')
        if np.sqrt(norm(pb1.ex-pb2.ex)**2\
                 + norm(pb1.ey-pb2.ey)**2\
                 + norm(pb1.ez-pb2.ez)**2) < 5e-3:
            print('Smaller that 5e-3, tolerated.')
        else:
            raise ValueError('Too large! Stopping.')
        
    # Check that there is no longitudinal separation
    try:
        assert(np.abs(np.dot(vbb_12, pb1.ez))<1e-4)
    except AssertionError:
        print(name_bb, 'The beams are longitudinally shifted')

    # Find separations
    sep_x.append(np.dot(vbb_12, pb1.ex))
    sep_y.append(np.dot(vbb_12, pb1.ey))
   
line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.lhcb1)

i_bb = 0
assert(bb_coupling==False)#Not implemented
for ee, eename in zip(line.elements, line.element_names):
    if isinstance(ee, pysixtrack.elements.BeamBeam4D):
        assert(eename==bb_names_b1[i_bb])
        
        ee.charge = bunch_intensity
        ee.sigma_x = np.sqrt(bb_sigmas_b2[11][i_bb]) 
        ee.sigma_y = np.sqrt(bb_sigmas_b2[33][i_bb])
        ee.beta_r = beta_r
        ee.x_bb = sep_x[i_bb]
        ee.y_bb = sep_y[i_bb]

        i_bb += 1
    if isinstance(ee, pysixtrack.elements.BeamBeam6D):
        assert(eename==bb_names_b1[i_bb])
        
        dpx = bb_xyz_b1[i_bb].tpx - bb_xyz_b2[i_bb].tpx
        dpy = bb_xyz_b1[i_bb].tpy - bb_xyz_b2[i_bb].tpy
        
        alpha, phi = find_alpha_and_phi(dpx, dpy)
        
        ee.phi = phi
        ee.alpha = alpha 
        ee.x_bb_co = sep_x[i_bb]
        ee.y_bb_co = sep_y[i_bb]
        ee.charge_slices = [bunch_intensity/n_slices]
        ee.zeta_slices = [0.0]
        ee.sigma_11 = bb_sigmas_b2[11][i_bb] 
        ee.sigma_12 = bb_sigmas_b2[12][i_bb]
        ee.sigma_13 = bb_sigmas_b2[13][i_bb]
        ee.sigma_14 = bb_sigmas_b2[14][i_bb] 
        ee.sigma_22 = bb_sigmas_b2[22][i_bb]
        ee.sigma_23 = bb_sigmas_b2[23][i_bb]
        ee.sigma_24 = bb_sigmas_b2[24][i_bb]
        ee.sigma_33 = bb_sigmas_b2[33][i_bb]
        ee.sigma_34 = bb_sigmas_b2[34][i_bb]  
        ee.sigma_44 = bb_sigmas_b2[44][i_bb]

        if not(bb_coupling):
            ee.sigma_13 = 0.
            ee.sigma_14 = 0.
            ee.sigma_23 = 0.
            ee.sigma_24 = 0.

        i_bb += 1

mad_ft = Madx()
mad_ft.call('mad/lhcwbb_fortracking.seq')
# without this the sequence does not work properly
mad_ft.use('lhcb1')

line_for_tracking, _ = pysixtrack.Line.from_madx_sequence(
        mad_ft.sequence['lhcb1'])


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


