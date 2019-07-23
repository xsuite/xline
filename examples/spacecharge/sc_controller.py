import numpy as np

import pysixtrack

class SC_controller(object):

    def __init__(self,type,intensity,eps_x,eps_y,dpp_rms,
                    bunchlength_rms=None,min_sigma_diff=1e-1):
        assert(type in ['SpaceChargeCoast', 'SpaceChargeBunched'])
        self._listSCnodes = []
        self._min_sigma_diff = min_sigma_diff
        self._type = type
        self.setBunchParameters(intensity,eps_x,eps_y,dpp_rms,bunchlength_rms)
    
    def installSCnodes(self,line,twiss,max_distance,centered=False):
        assert (len(twiss['s']) == len(line))
        gamma = twiss.summary.gamma
        beta = np.sqrt(1.-1./gamma**2)
        self._circumference = twiss.summary.length

        circumference = self._circumference
        intensity = self._intensity
        eps_x = self._eps_x
        eps_y = self._eps_y
        dpp_rms = self._dpp_rms
        bunchlength_rms = self._bunchlength_rms

        new_line = pysixtrack.Line(elements=[], element_names=[])
        ss_lastSCkick = 0
        for ii, (ee, nn)  in enumerate(zip(line.elements, line.element_names)):
            ss = twiss['s'][ii]
            new_line.append_element(ee, nn)

            # check what the next length would be
            jj = ii
            ss_next = ss
            while (jj<len(line) and ss_next==ss):
                ss_next = twiss['s'][jj]
                jj += 1
            predicted_length = ss_next - ss_lastSCkick
            if predicted_length>max_distance or ii==(len(line)-1):
                beta_x = twiss['betx'][ii]
                beta_y = twiss['bety'][ii]
                D_x = twiss['dx'][ii]
                D_y = twiss['dy'][ii]
                CO_x = twiss['x'][ii]
                CO_y = twiss['y'][ii]

                sigma_x = self._sigma(beta_x, eps_x, D_x, dpp_rms)
                sigma_y = self._sigma(beta_y, eps_y, D_y, dpp_rms)
                length = ss-ss_lastSCkick
                if length > max_distance:
                    print("Warning: length of space charge node bigger than max_distance")
                if self._type == 'SpaceChargeCoast':
                    scNode = pysixtrack.elements.SpaceChargeCoast(
                        line_density=intensity/circumference,
                        sigma_x=sigma_x,
                        sigma_y=sigma_y,
                        length=length,
                        min_sigma_diff=self._min_sigma_diff,
                        Delta_x=CO_x,
                        Delta_y=CO_y,
                        enabled=True)
                if self._type == 'SpaceChargeBunched':
                    scNode = pysixtrack.elements.SpaceChargeBunched(
                        number_of_particles=intensity,
                        bunchlength_rms=bunchlength_rms,
                        sigma_x=sigma_x,
                        sigma_y=sigma_y,
                        length=length,
                        min_sigma_diff=self._min_sigma_diff,
                        Delta_x=CO_x,
                        Delta_y=CO_y,
                        enabled=True)
                new_line.append_element(scNode, 'SC%i'%(len(self._listSCnodes)))
                ss_lastSCkick = ss

                optics = {}
                optics['beta_x'] = beta_x
                optics['beta_y'] = beta_y
                optics['D_x'] = D_x
                optics['D_y'] = D_y
                optics['s'] = ss
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

        return new_line

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


