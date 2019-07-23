import pickle

import numpy as np
import matplotlib.pylab as plt

from pyoptics import madlang,optics
import pyoptics.tfsdata as tfs
import pysixtrack
from pysixtrack.particles import Particles
import pysixtrack.helpers as hp 

from sc_controller import SC_controller

#see sps/madx/a001_track_thin.madx
mad=madlang.open('madx/SPS_Q20_thin.seq')
mad.acta_31637.volt=4.5
mad.acta_31637.lag=0.5

dict_elements = {}
for ee in dir(pysixtrack.elements):
    if ee.startswith('_'):
        continue
    if ee == 'np':
        continue
    dict_elements[ee] = getattr(pysixtrack.elements, ee)


elems,rest,iconv=mad.sps.expand_struct(dict_elements)
sps=pysixtrack.Line(elements= [e[2] for e in elems])
spstwiss=tfs.open('madx/twiss_SPS_Q20_thin.tfs')

# remove beg. and end. markers in twiss table (as they are not in sixtrack lattice)
for k in spstwiss:
    try:
        spstwiss[k] = spstwiss[k][1:-1]
    except TypeError:
        pass
assert (len(spstwiss['s']) == len(elems))

gamma = spstwiss['param']['gamma']
beta = np.sqrt(1.-1./gamma**2)
betagamma = beta*gamma

p0c = 25.92e9
intensity=2e11 #* spstwiss['param']['length'] # for a coasting beam with a line density of 2e11/m
eps_x=2e-6/betagamma
eps_y=2e-6/betagamma
dpp_rms=1.5e-3
bunchlength_rms = 0.22

max_distance = 35. #6.9 #25.
# my_SC_controller = SC_controller('SpaceChargeCoast',intensity * spstwiss['param']['length'],eps_x,eps_y,dpp_rms)
my_SC_controller = SC_controller('SpaceChargeBunched',intensity,eps_x,eps_y,dpp_rms,bunchlength_rms)
new_elems = my_SC_controller.installSCnodes(elems,spstwiss,max_distance=max_distance,centered=False) #25.

# prepare a particle on the closed orbit
p=Particles(p0c=p0c)
ring = hp.Ring(new_elems, p0c=p0c)

my_SC_controller.disableSCnodes()
# print("\n  starting closed orbit search ... ")
closed_orbit = ring.find_closed_orbit(guess=[spstwiss['x'][0], spstwiss['px'][0], 
    spstwiss['y'][0], spstwiss['py'][0], 0., 0.], method='get_guess')
my_SC_controller.enableSCnodes()

with open('particle_on_CO.pkl', 'wb') as fid:
    closed_orbit[0]._m = None # to be sorted out 
    pickle.dump(closed_orbit[0], fid)
with open('line.pkl', 'wb') as fid:
    pickle.dump(new_elems, fid)
with open('SCcontroller.pkl', 'wb') as fid:
    pickle.dump(my_SC_controller, fid)



import matplotlib.patches as patches
if 1:
    plt.close('all')

    f, ax = plt.subplots()
    sc_lengths = my_SC_controller.getSCnodesLengths()
    ax.hist(sc_lengths, bins=np.linspace(0,max(sc_lengths)+0.1,100))
    ax.set_xlabel('length of SC kick (m)')
    ax.set_ylabel('counts')
    ax.set_xlim(left=0)
    plt.show()

    f, ax = plt.subplots(figsize=(14,5))
    ax.plot(spstwiss['s'], spstwiss['betx'], 'b', label='x', lw=2)
    ax.plot(spstwiss['s'], spstwiss['bety'], 'g', label='y', lw=2)
    for s in my_SC_controller.getSCnodesOpticsParameter('s'): ax.axvline(s, linewidth=1, color='r', linestyle='--')
    ax.set_xlim(0,1100)
    ax.set_ylim(0,120)
    ax.set_xlabel('s (m)')
    ax.set_ylabel('beta functions (m)')
    ax.legend(loc=3)
    plt.show()


