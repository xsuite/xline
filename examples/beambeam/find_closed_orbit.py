import sixtracktools
import pysixtrack
import pysixtrack.helpers as hp
from pysixtrack.particles import Particles

import scipy.optimize as so

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

        pcl=Particles(p0c=self.p0c)

        pcl.x = coord[0]
        pcl.px = coord[1]
        pcl.y = coord[2]
        pcl.py = coord[3]
        pcl.zeta = coord[4]
        pcl.delta = coord[5]

        self.track(pcl)

        coord_out = np.array([pcl.x, pcl.px, pcl.y, pcl.py, pcl.sigma, pcl.delta])

        return coord_out

    def _CO_error(self, coord):
        return np.sum((self.one_turn_map(coord)-coord)**2)

    def find_closed_orbit(self, guess=np.array([0.,0.,0.,0.,0.,0.])):
        res = so.minimize(self._CO_error, guess, tol=1e-20, method='Nelder-Mead')

        pcl=Particles(p0c=self.p0c)

        pcl.x = res.x[0]
        pcl.px = res.x[1]
        pcl.y = res.x[2]
        pcl.py = res.x[3]
        pcl.zeta = res.x[4]
        pcl.delta = res.x[5]

        CO = self.track_elem_by_elem(pcl)

        return CO

ring = Ring(line, p0c=p0c_eV)

closed_orbit = ring.find_closed_orbit()

x_CO = np.array([pp.x for pp in closed_orbit])
y_CO = np.array([pp.y for pp in closed_orbit])

# Compare closed orbit against sixtrack
# Load sixtrack tracking data
sixdump_all = sixtracktools.SixDump101('res/dump3.dat')
# Assume first particle to be on the closed orbit
Nele_st = len(iconv)
sixdump_CO = sixdump_all[::2][:Nele_st]
x_CO_at_st_ele = x_CO[iconv]
y_CO_at_st_ele = y_CO[iconv]
print('Max C.O. discrepancy in x %.2e m'%np.max(np.abs(x_CO_at_st_ele-sixdump_CO.x)))
print('Max C.O. discrepancy in y %.2e m'%np.max(np.abs(y_CO_at_st_ele-sixdump_CO.y)))