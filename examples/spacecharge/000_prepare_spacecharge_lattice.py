from pyoptics import madlang,optics
import pyoptics.tfsdata as tfs
import pysixtrack
from pysixtrack.particles import Particles
import pysixtrack.helpers as hp 
import numpy as np
import matplotlib.pylab as plt
import pickle

class SC_controller(object):

    def __init__(self,type,intensity,eps_x,eps_y,dpp_rms,
                    bunchlength_rms=None,min_sigma_diff=1e-1):
        assert(type in ['SpaceChargeCoast', 'SpaceChargeBunched'])
        self._listSCnodes = []
        self._min_sigma_diff = min_sigma_diff
        self._type = type
        self.setBunchParameters(intensity,eps_x,eps_y,dpp_rms,bunchlength_rms)
    
    def installSCnodes(self,lattice,twiss,min_distance,centered=False):
        assert (len(twiss['s']) == len(lattice))
        gamma = twiss['param']['gamma']
        beta = np.sqrt(1.-1./gamma**2)
        self._circumference = twiss['param']['length']

        circumference = self._circumference
        intensity = self._intensity
        eps_x = self._eps_x
        eps_y = self._eps_y
        dpp_rms = self._dpp_rms
        bunchlength_rms = self._bunchlength_rms

        new_lattice = []
        s_lastSCkick = 0
        for i, e in enumerate(lattice):
            s = twiss['s'][i]
            new_lattice.append(e)
            if s-s_lastSCkick >= min_distance or i==(len(elems)-1):
                beta_x = twiss['betx'][i]
                beta_y = twiss['bety'][i]
                D_x = twiss['dx'][i]
                D_y = twiss['dy'][i]
                CO_x = twiss['x'][i]
                CO_y = twiss['y'][i]

                sigma_x = self._sigma(beta_x, eps_x, D_x, dpp_rms)
                sigma_y = self._sigma(beta_y, eps_y, D_y, dpp_rms)
                length = s-s_lastSCkick
                if self._type == 'SpaceChargeCoast':
                    scNode = pysixtrack.SpaceChargeCoast(
                        line_density=intensity/circumference,
                        sigma_x=sigma_x,
                        sigma_y=sigma_y,
                        length=length,
                        min_sigma_diff=self._min_sigma_diff,
                        Delta_x=CO_x,
                        Delta_y=CO_y,
                        enabled=True)
                if self._type == 'SpaceChargeBunched':
                    scNode = pysixtrack.SpaceChargeBunched(
                        number_of_particles=intensity,
                        bunchlength_rms=bunchlength_rms,
                        sigma_x=sigma_x,
                        sigma_y=sigma_y,
                        length=length,
                        min_sigma_diff=self._min_sigma_diff,
                        Delta_x=CO_x,
                        Delta_y=CO_y,
                        enabled=True)
                new_lattice.append(('SC%i'%(len(self._listSCnodes)), self._type, scNode))
                s_lastSCkick = s

                optics = {}
                optics['beta_x'] = beta_x
                optics['beta_y'] = beta_y
                optics['D_x'] = D_x
                optics['D_y'] = D_y
                optics['s'] = s
                self._addNode(optics, scNode)
        
        if centered:
            sc_nodes = self.getSCnodes()[0]
            for i in range(len(sc_nodes)-1):
                new_length = (sc_nodes[i].length + sc_nodes[i+1].length)/2.
                sc_nodes[i].length = new_length
            sc_nodes[-1].length += circumference - sum(self.getSCnodesLengths())

        lengthsSCndodes = self.getSCnodesLengths()
        print("\n  number of installed space charge kicks: %d"%len(self._listSCnodes))
        print("\n  maximum length of space charge kicks: %1.2f m"%(max(lengthsSCndodes)))
        print("\n  average length of space charge kicks: %1.2f m"%(np.mean(lengthsSCndodes)))
        print("\n  integrated length of space charge kicks: %1.4f m"%(sum(lengthsSCndodes)))
        
        # verify that integrated length of SC kicks corresponds to circumference
        assert(self._circumference == sum(lengthsSCndodes))

        return new_lattice

    def getSCnodesOpticsParameter(self,p):
        return [optics[p] for optics in self.getSCnodes()[1]]

    # def setSCnodesOpticsParameter(self,p,newValues):
    #     assert(len(newValues) == len(self._listSCnodes))
    #     for i, optics in enumerate(self.getSCnodes()[1]):
    #         optics[p] = newValues[i]
    
    def setClosedOrbit(self):
        print('... to be implemented ...')

    def _sigma(self, beta, eps, D, dpp_rms):
        return np.sqrt(beta*eps + D**2*dpp_rms**2)
    
    def _addNode(self, node, optics):
        self._listSCnodes.append((optics, node))
    
    def getSCnodes(self):
        # nodes, optics = list(zip(*self.listSCnodes))
        nodes = [sc[0] for sc in self._listSCnodes]
        optics = [sc[1] for sc in self._listSCnodes]
        return nodes, optics

    def getSCnodesLengths(self):
        nodes, _ = self.getSCnodes()
        return [n.length for n in nodes]
    
    def disableSCnodes(self):
        nodes, _ = self.getSCnodes()
        for n in nodes: n.enabled = False
    
    def enableSCnodes(self):
        nodes, _ = self.getSCnodes()
        for n in nodes: n.enabled = True
    
    def setBunchParameters(self,intensity,eps_x,eps_y,dpp_rms,bunchlength_rms=None):
        self._intensity = intensity
        self._eps_x = eps_x
        self._eps_y = eps_y
        self._dpp_rms = dpp_rms
        self._bunchlength_rms = bunchlength_rms
        if self._listSCnodes:
            for node, optics in zip(self.getSCnodes()):
                beta_x = optics['beta_x']
                beta_y = optics['beta_y']
                D_x = optics['D_x']
                D_y = optics['D_y']

                node.sigma_x = self._sigma(beta_x, eps_x, D_x, dpp_rms)
                node.sigma_y = self._sigma(beta_y, eps_y, D_y, dpp_rms)
                if self._type == 'SpaceChargeCoast':
                    node.line_density = intensity/self._circumference
                if self._type == 'SpaceChargeBunched':
                    node.number_of_particles = intensity
                    node.bunchlength_rms = bunchlength_rms



#see sps/madx/a001_track_thin.madx
mad=madlang.open('madx/SPS_Q20_thin.seq')
mad.acta_31637.volt=4.5
mad.acta_31637.lag=0.5


elems,rest,iconv=mad.sps.expand_struct(pysixtrack.element_types)
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

min_distance = 5. #6.9 #25.
# my_SC_controller = SC_controller('SpaceChargeCoast',intensity * spstwiss['param']['length'],eps_x,eps_y,dpp_rms)
my_SC_controller = SC_controller('SpaceChargeBunched',intensity,eps_x,eps_y,dpp_rms,bunchlength_rms)
new_elems = my_SC_controller.installSCnodes(elems,spstwiss,min_distance=min_distance,centered=False) #25.

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


plot = 1 #False
import matplotlib.patches as patches
if plot:
    plt.close('all')

    f, ax = plt.subplots()
    ax.hist(my_SC_controller.getSCnodesLengths(), 100)
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


