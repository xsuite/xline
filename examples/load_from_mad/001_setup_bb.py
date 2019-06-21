import numpy as np
import pickle

_sigma_names = [11,12,13,14,22,23,24,33,34,44] 

class MadPoint(object):
    def __init__(self,name,mad, add_CO=True):
        self.name=name
        twiss=mad.table.twiss
        survey=mad.table.survey
        idx=np.where(survey.name==name)[0][0]
        self.tx=twiss.x[idx]
        self.ty=twiss.y[idx]
        self.tpx=twiss.px[idx]
        self.tpy=twiss.py[idx]
        self.sx=survey.x[idx]
        self.sy=survey.y[idx]
        self.sz=survey.z[idx]
        theta=survey.theta[idx]
        phi=survey.phi[idx]
        psi=survey.psi[idx]
        thetam=np.array([[np.cos(theta) ,           0,np.sin(theta)],
             [          0,           1,         0],
             [-np.sin(theta),           0,np.cos(theta)]])
        phim=np.array([[          1,          0,          0],
            [          0,np.cos(phi)   ,   np.sin(phi)],
            [          0,-np.sin(phi)  ,   np.cos(phi)]])
        psim=np.array([[   np.cos(psi),  -np.sin(psi),          0],
            [   np.sin(psi),   np.cos(psi),          0],
            [          0,          0,          1]])
        wm=np.dot(thetam,np.dot(phim,psim))
        self.ex=np.dot(wm,np.array([1,0,0]))
        self.ey=np.dot(wm,np.array([0,1,0]))
        self.ez=np.dot(wm,np.array([0,0,1]))
        self.sp=np.array([self.sx,self.sy,self.sz])
        if add_CO:
            self.p=self.sp+ self.ex * self.tx + self.ey * self.ty
        else:
            self.p=self.sp

    def dist(self,other):
        return np.sqrt(np.sum((self.p-other.p)**2))

    def distxy(self,other):
        dd=self.p-other.p
        return np.dot(dd,self.ex),np.dot(dd,self.ey)

def find_alpha_and_phi(dpx, dpy):
    
    phi = np.sqrt(dpx**2+dpy**2)/2.
    if phi < 1e-20:
        alpha = 0.
    elif np.abs(dpx) >= np.abs(dpy):
        alpha = np.arctan(dpy/dpx)
        if dpx<0: phi = -phi
    else:
        alpha = np.sign(dpy)*(np.pi/2-np.abs(np.arctan(dpx/dpy)))
        if dpy<0: phi = -phi
    
    return alpha, phi


def get_bb_names_xyz_points_sigma_matrices(mad, seq_name):
    mad.use(sequence=seq_name);
    mad.twiss()
    mad.survey()
    
    seq = mad.sequence[seq_name]
   
    bb_names = []
    bb_xyz_points = []
    bb_sigmas = {kk:[] for kk in _sigma_names}
    
    for ee in seq.elements:
        if ee.base_type.name == 'beambeam':
            eename = ee.name
            bb_names.append(eename)
            bb_xyz_points.append(MadPoint(eename+':1', mad))

            i_twiss = np.where(mad.table.twiss.name==(eename+':1'))[0][0]
            
            for sn in _sigma_names:
                bb_sigmas[sn].append(
                        getattr(mad.table.twiss, 'sig%d'%sn)[i_twiss])

    return bb_names, bb_xyz_points, bb_sigmas

def norm(v):
    return np.sqrt(np.sum(v**2))


from cpymad.madx import Madx
import pysixtrack

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

# IP locations
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

# Check number of slice
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
            print('Smaller that 5e-3, tollerated.')
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
   
    if name_bb == 'bb_par.l5b1_21':
        prrrr



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
        ee.charge_slices = [0.0]
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

dct = line.to_dict(keepextra=True)
#dct['elements'] = dct['elements'][:2]
with open('line_from_mad.pkl', 'wb') as fid:
    pickle.dump(dct, fid)

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


