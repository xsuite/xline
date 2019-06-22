import numpy as np
import pysixtrack
from pysixtrack import MadPoint

_sigma_names = [11,12,13,14,22,23,24,33,34,44] 

def norm(v):
    return np.sqrt(np.sum(v**2))

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


def get_bb_names_madpoints_sigmas(mad, seq_name):
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


def shift_strong_beam_based_on_close_ip(points_weak, points_strong,
        IPs_survey_weak, IPs_survey_strong):

    for i_bb, _ in enumerate(points_weak):
        
        pbw = points_weak[i_bb]
        pbs = points_strong[i_bb]
        
        # Find closest IP
        d_ip = 1e6
        use_ip = 0
        for ip in IPs_survey_weak.keys():
            dd = norm(pbw.p - IPs_survey_weak[ip].p)
            if dd < d_ip:
                use_ip = ip
                d_ip = dd
    
        # Shift Bs
        shift_ws = IPs_survey_strong[use_ip].p - IPs_survey_weak[use_ip].p
        pbs.p -= shift_ws
    
def find_bb_separations(points_weak, points_strong,names=None):

    if names is None:
        names = ['bb_%d'%ii for ii in range(len(points_weak))]

    sep_x = []
    sep_y = []
    for i_bb, name_bb in enumerate(names):
        
        pbw = points_weak[i_bb]
        pbs = points_strong[i_bb]
    
        # Find vws
        vbb_ws = points_strong[i_bb].p - points_weak[i_bb].p
    
        # Check that the two reference system are parallel
        try:
            assert(norm(pbw.ex-pbs.ex)<1e-10) #1e-4 is a reasonable limit
            assert(norm(pbw.ey-pbs.ey)<1e-10) #1e-4 is a reasonable limit
            assert(norm(pbw.ez-pbs.ez)<1e-10) #1e-4 is a reasonable limit
        except AssertionError:
            print(name_bb, 'Reference systems are not parallel')
            if np.sqrt(norm(pbw.ex-pbs.ex)**2\
                     + norm(pbw.ey-pbs.ey)**2\
                     + norm(pbw.ez-pbs.ez)**2) < 5e-3:
                print('Smaller that 5e-3, tolerated.')
            else:
                raise ValueError('Too large! Stopping.')
            
        # Check that there is no longitudinal separation
        try:
            assert(np.abs(np.dot(vbb_ws, pbw.ez))<1e-4)
        except AssertionError:
            print(name_bb, 'The beams are longitudinally shifted')
        
        # Find separations
        sep_x.append(np.dot(vbb_ws, pbw.ex))
        sep_y.append(np.dot(vbb_ws, pbw.ey))
   
    return sep_x, sep_y

def setup_beam_beam_in_line(line, bb_names, bb_sigmas_strong, 
    bb_points_weak, bb_points_strong, beta_r_strong,
    bunch_intensity_strong, n_slices_6D, bb_coupling):

    bb_sep_x, bb_sep_y = find_bb_separations(
        points_weak=bb_points_weak, points_strong=bb_points_strong,
        names = bb_names) 

    i_bb = 0
    assert(bb_coupling==False)#Not implemented
    for ee, eename in zip(line.elements, line.element_names):
        if isinstance(ee, pysixtrack.elements.BeamBeam4D):
            assert(eename==bb_names[i_bb])
            
            ee.charge = bunch_intensity_strong
            ee.sigma_x = np.sqrt(bb_sigmas_strong[11][i_bb]) 
            ee.sigma_y = np.sqrt(bb_sigmas_strong[33][i_bb])
            ee.beta_r = beta_r
            ee.x_bb = sep_x[i_bb]
            ee.y_bb = sep_y[i_bb]
    
            i_bb += 1
        if isinstance(ee, pysixtrack.elements.BeamBeam6D):
            assert(eename==bb_names[i_bb])
            
            dpx = bb_points_weak[i_bb].tpx - bb_points_strong[i_bb].tpx
            dpy = bb_points_weak[i_bb].tpy - bb_points_strong[i_bb].tpy
            
            alpha, phi = find_alpha_and_phi(dpx, dpy)
            
            ee.phi = phi
            ee.alpha = alpha 
            ee.x_bb_co = sep_x[i_bb]
            ee.y_bb_co = sep_y[i_bb]
            ee.charge_slices = [bunch_intensity/n_slices]
            ee.zeta_slices = [0.0]
            ee.sigma_11 = bb_sigmas_strong[11][i_bb] 
            ee.sigma_12 = bb_sigmas_strong[12][i_bb]
            ee.sigma_13 = bb_sigmas_strong[13][i_bb]
            ee.sigma_14 = bb_sigmas_strong[14][i_bb] 
            ee.sigma_22 = bb_sigmas_strong[22][i_bb]
            ee.sigma_23 = bb_sigmas_strong[23][i_bb]
            ee.sigma_24 = bb_sigmas_strong[24][i_bb]
            ee.sigma_33 = bb_sigmas_strong[33][i_bb]
            ee.sigma_34 = bb_sigmas_strong[34][i_bb]  
            ee.sigma_44 = bb_sigmas_strong[44][i_bb]
    
            if not(bb_coupling):
                ee.sigma_13 = 0.
                ee.sigma_14 = 0.
                ee.sigma_23 = 0.
                ee.sigma_24 = 0.
    
            i_bb += 1


