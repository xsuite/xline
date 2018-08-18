import sixtracktools
import pysixtrack
import pysixtrack.helpers as hp
from pysixtrack.particles import Particles

import numpy as np

six = sixtracktools.SixTrackInput('.')
p0c_eV = six.initialconditions[-3]*1e6

line, rest, iconv = six.expand_struct(convert=pysixtrack.element_types)

ind_BB4D, namelistBB4D, listBB4D = hp.get_elems_of_type(line, 'BeamBeam4D')

for bb in listBB4D:
    bb.enabled = False



class Ring(pysixtrack.Line):
    def __init__(self, line, p0c=1e9):
        self.p0c = p0c
        super().__init__(elements=[elem for label,elem_type,elem in line])

    def one_turn_map(self, coord):

        pcl=Particles(p0c=self.p0c_eV)

        pcl.x = coord[0]
        pcl.px = coord[1]
        pcl.y = coord[2]
        pcl.py = coord[3]
        pcl.sigma = coord[4]
        pcl.delta = coord[5]

        self.track(pcl)

        coord_out = np.array([pcl.x, pcl.px, pcl.y, pcl.py, pcl.sigma, pcl.delta])

        return coord_out

ring = Ring(line, p0c=p0c_eV)

prrrr

# Build 1-turn map function
def one_turn_map(coord_in, flag_ebe=False):

    # coord = np.array([x, px, y, py, sigma, delta])

    coord = coord_in
    
    npart = 1
    
    bunch=sixtracklib.CParticles(npart=npart, 
            p0c=p0c_eV,
            beta0 = beta0,
            gamma0 = gamma0)
    bunch.x+=coord[0]
    bunch.px+=coord[1]
    bunch.y+=coord[2]
    bunch.py+=coord[3]
    bunch.sigma+=coord[4]
    bunch.set_delta(coord[5])

    particles,ebe,tbt = machine.track_cl(bunch,nturns=1,elembyelem={True:True, False:None}[flag_ebe], turnbyturn=True)

    coord =  np.array([tbt.x[1][0], tbt.px[1][0], tbt.y[1][0], tbt.py[1][0], 
                    tbt.sigma[1][0], tbt.delta[1][0]])

    if flag_ebe:
        return coord, ebe
    else:
        return coord

# Define function for optimization
tominimize = lambda coord: np.sum((one_turn_map(coord)-coord)**2)

# Find fixed point
res = so.minimize(tominimize, np.array([0.,0.,0.,0.,0.,0.]), tol=1e-20, method='Nelder-Mead')
# Get close orbit around the machine
temp, ebe_CO = one_turn_map(res.x, flag_ebe=True)


# sixdump_all = sixtracktools.SixDump101('res/dump3.dat')



