import numpy as np

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


def shift_strong_beam_based_on_closest_ip(points_bw, points_bs,
        IPs_survey_bw, IPs_survey_bs):

    for i_bb, _ in enumerate(points_bw):
        
        pbw = points_bw[i_bb]
        pbs = points_bs[i_bb]
        
        # Find closest IP
        d_ip = 1e6
        use_ip = 0
        for ip in ip_names:
            dd = norm(pbw.p - IP_xyz_bw[ip].p)
            if dd < d_ip:
                use_ip = ip
                d_ip = dd
    
        # Shift Bs
        shift_ws = IP_survey_bs[use_ip].p - IP_survey_bw[use_ip].p
        pbs.p -= shift_ws
    
           
 


