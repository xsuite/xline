import sixtracktools
import pysixtrack
import pysixtrack.helpers as hp


import numpy as np

six = sixtracktools.SixTrackInput('.')
p0c_eV = six.initialconditions[-3]*1e6

line, rest, iconv = six.expand_struct(convert=pysixtrack.element_types)

# Disable BB elements
ind_BB4D, namelistBB4D, listBB4D = hp.get_elems_of_type(line, 'BeamBeam4D')
for bb in listBB4D: bb.enabled = False

#Find closed orbit
ring = hp.Ring(line, p0c=p0c_eV)
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