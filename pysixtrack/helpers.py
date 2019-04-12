import numpy as np
import pysixtrack

import scipy.optimize as so

from pysixtrack.particles import Particles


def get_elems_of_type(line, elemname):

    pyst_ele_type_list = [ele[1] for ele in line]
    indlist = np.where([ss == elemname for ss in pyst_ele_type_list])[0]
    ele_list = [line[ind][2] for ind in indlist]
    name_list = [line[ind][0] for ind in indlist]

    return indlist, name_list, ele_list


class Ring(pysixtrack.Line):
    def __init__(self, line, p0c=1e9):
        self.p0c = p0c
        super().__init__(elements=[elem for label, elem_type, elem in line])

    def one_turn_map(self, coord):

        pcl = Particles(p0c=self.p0c)

        pcl.x = coord[0]
        pcl.px = coord[1]
        pcl.y = coord[2]
        pcl.py = coord[3]
        pcl.zeta = coord[4]
        pcl.delta = coord[5]

        self.track(pcl)

        coord_out = np.array(
            [pcl.x, pcl.px, pcl.y, pcl.py, pcl.sigma, pcl.delta])

        return coord_out

    def _CO_error(self, coord):
        return np.sum((self.one_turn_map(coord) - coord)**2)

    def find_closed_orbit(
            self,
            guess=[
                0.,
                0.,
                0.,
                0.,
                0.,
                0.],
            method='Nelder-Mead'):

        if method == 'get_guess':
            res = type('', (), {})()
            res.x = guess
        else:
            res = so.minimize(self._CO_error, np.array(
                guess), tol=1e-20, method=method)

        pcl = Particles(p0c=self.p0c)

        pcl.x = res.x[0]
        pcl.px = res.x[1]
        pcl.y = res.x[2]
        pcl.py = res.x[3]
        pcl.zeta = res.x[4]
        pcl.delta = res.x[5]

        CO = self.track_elem_by_elem(pcl)

        return CO
