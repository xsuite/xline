import pickle

import numpy as np
import matplotlib.pylab as plt

from cpymad.madx import Madx

import pysixtrack
from pysixtrack.particles import Particles

from sc_controller import SC_controller


mad = Madx()
mad.options.echo=False
mad.options.info=False
mad.warn=False 

mad.call('madx/SPS_Q20_thin.seq')
mad.use('sps')

mad.globals.klsfa = 0.
mad.globals.klsfb = 0.
mad.globals.klsfc = 0.
mad.globals.klsda = 0.
mad.globals.klsdb = 0.

twtable = mad.twiss()

mad.elements['acta.31637'].volt = 0*4.5
mad.elements['acta.31637'].lag = 0.5


line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.sps)

# Checking our assumptions
assert(len(twtable.name) == len(line.element_names))
for tnn, lnn in zip(twtable.name, line.element_names):
    assert(tnn.split(':')[0] == lnn)

gamma = twtable.summary.gamma 
beta = np.sqrt(1.-1./gamma**2)
betagamma = beta*gamma

p0c = 25.92e9
intensity=0*2e11 #* spstwiss['param']['length'] # for a coasting beam with a line density of 2e11/m
eps_x=2e-6/betagamma
eps_y=2e-6/betagamma
dpp_rms=1.5e-3
bunchlength_rms = 0.22

max_distance = 35. #6.9 #25.
# my_SC_controller = SC_controller('SpaceChargeCoast',intensity * spstwiss['param']['length'],eps_x,eps_y,dpp_rms)
my_SC_controller = SC_controller('SpaceChargeBunched',intensity,eps_x,eps_y,dpp_rms,bunchlength_rms)
new_elems = my_SC_controller.installSCnodes(line,twtable,max_distance=max_distance,centered=False) #25.



# prepare a particle on the closed orbit
# p=Particles(p0c=p0c)

my_SC_controller.disableSCnodes()
# print("\n  starting closed orbit search ... ")
part_on_CO = line.find_closed_orbit(guess=[twtable['x'][0], twtable['px'][0], 
    twtable['y'][0], twtable['py'][0], 0., 0.], p0c=p0c, method='get_guess')
closed_orbit = line.track_elem_by_elem(part_on_CO)
#my_SC_controller.enableSCnodes()

with open('particle_on_CO.pkl', 'wb') as fid:
    closed_orbit[0]._m = None # to be sorted out 
    pickle.dump(closed_orbit[0], fid)
with open('line.pkl', 'wb') as fid:
    pickle.dump(new_elems, fid)
with open('SCcontroller.pkl', 'wb') as fid:
    pickle.dump(my_SC_controller, fid)

with open('twiss_at_start.pkl', 'wb') as fid:
    pickle.dump({
        'betx': twtable.betx[0],
        'bety': twtable.bety[0]}, fid)

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
    ax.plot(twtable.s, twtable.betx, 'b', label='x', lw=2)
    ax.plot(twtable.s, twtable.bety, 'g', label='x', lw=2)
    for s in my_SC_controller.getSCnodesOpticsParameter('s'): ax.axvline(s, linewidth=1, color='r', linestyle='--')
    ax.set_xlim(0,1100)
    ax.set_ylim(0,120)
    ax.set_xlabel('s (m)')
    ax.set_ylabel('beta functions (m)')
    ax.legend(loc=3)
    plt.show()


