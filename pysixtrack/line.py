import numpy as np

from .base_classes import Element
from . import elements
from .particles import Particles

from .loader_sixtrack import _expand_struct
from .loader_mad import _from_madx_sequence


class Line(Element):
    _description = [
        ('elements', '', 'List of elements', ()),
        ('element_names', '', 'List of element names', ())
    ]
    _extra = []

    def to_dict(self, keepextra=False):
        out = {}
        out['elements'] = [el.to_dict(keepextra) for el in self.elements]
        out['element_names'] = self.element_names[:]
        return out

    @classmethod
    def from_dict(cls,dct,keepextra=True):
        self=cls(elements=[], element_names=[])
        for el in dct['elements']:
            eltype=getattr(elements,el['__class__'])
            newel=eltype.from_dict(el)
            self.elements.append(newel)
        self.element_names=dct['element_names']
        return self

    def append_line(self, line):
        # Append the elements
        if type(line) is Line:
            # got a pysixtrack line
            self.elements += line.elements
        else:
            # got a different type of line (e.g. pybplep)
            for ee in line.elements:
                type_name = ee.__class__.__name__
                newele = element_types[type_name](**ee._asdict())
                self.elements.append(newele)

        # Append the names
        self.element_names += line.element_names

        assert(len(self.elements) == len(self.element_names))

        return self

    def track(self, p):
        for el in self.elements:
            el.track(p)

    def track_elem_by_elem(self, p):
        out = []
        for el in self.elements:
            out.append(p.copy())
            el.track(p)
        return out

    def find_closed_orbit(self, p0c, guess=[0., 0., 0., 0., 0., 0.],
                          method='Nelder-Mead'):

        def _one_turn_map(coord):
            pcl = Particles(p0c=p0c)
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

        def _CO_error(coord):
            return np.sum((_one_turn_map(coord) - coord)**2)

        if method == 'get_guess':
            res = type('', (), {})()
            res.x = guess
        else:
            import scipy.optimize as so
            res = so.minimize(_CO_error, np.array(
                guess), tol=1e-20, method=method)

        pcl = Particles(p0c=p0c)

        pcl.x = res.x[0]
        pcl.px = res.x[1]
        pcl.y = res.x[2]
        pcl.py = res.x[3]
        pcl.zeta = res.x[4]
        pcl.delta = res.x[5]

        return pcl

    def get_length(self):
        thick_element_types = (elements.Drift, elements.DriftExact)

        ll = 0
        for ee in self.elements:
            if isinstance(ee, thick_element_types):
                ll+=ee.length
        return ll
    
    def remove_zero_length_drifts(self):
        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            
            if isinstance(ee, (elements.Drift, elements.DriftExact)):
                if ee.length==0.:
                    continue

            newline.elements.append(ee)
            newline.element_names.append(nn)

        return newline

    def get_elements_of_type(self, types):
        if not hasattr(types, '__iter__'):
            type_list = [types]
        else:
            type_list = types

        names = []
        elements = []
        for ee, nn in zip(self.elements, self.element_names):
            for tt in type_list:
                if isinstance(ee, tt):
                    names.append(nn)
                    elements.append(ee)

        return elements, names



    def enable_beambeam(self):

        for ee in self.elements:
            if isinstance(ee, (elements.BeamBeam4D, elements.BeamBeam6D)):
                ee.enabled = True

    def disable_beambeam(self):

        for ee in self.elements:
            if isinstance(ee, (elements.BeamBeam4D, elements.BeamBeam6D)):
                ee.enabled = False

    def beambeam_store_closed_orbit_and_dipolar_kicks(self, particle_on_CO,
                                                      separation_given_wrt_closed_orbit_4D=True,
                                                      separation_given_wrt_closed_orbit_6D=True):

        self.disable_beambeam()
        closed_orbit = self.track_elem_by_elem(particle_on_CO)

        self.enable_beambeam()

        for ie, ee in enumerate(self.elements):
            # to transfer to beambeam.py

            if ee.__class__.__name__ == 'BeamBeam4D':
                if separation_given_wrt_closed_orbit_4D:
                    ee.x_bb += closed_orbit[ie].x
                    ee.y_bb += closed_orbit[ie].y

                # Evaluate dipolar kick
                ptemp = closed_orbit[ie].copy()
                ptempin = ptemp.copy()

                ee.track(ptemp)

                ee.d_px = ptemp.px - ptempin.px
                ee.d_py = ptemp.py - ptempin.py

            elif ee.__class__.__name__ == 'BeamBeam6D':
                if not separation_given_wrt_closed_orbit_6D:
                    raise ValueError('Not implemented!')

                # Store closed orbit
                ee.x_co = closed_orbit[ie].x
                ee.px_co = closed_orbit[ie].px
                ee.y_co = closed_orbit[ie].y
                ee.py_co = closed_orbit[ie].py
                ee.zeta_co = closed_orbit[ie].zeta
                ee.delta_co = closed_orbit[ie].delta

                # Evaluate 6d kick on closed orbit
                ptemp = closed_orbit[ie].copy()
                ptempin = ptemp.copy()

                ee.track(ptemp)

                ee.d_x = ptemp.x - ptempin.x
                ee.d_px = ptemp.px - ptempin.px
                ee.d_y = ptemp.y - ptempin.y
                ee.d_py = ptemp.py - ptempin.py
                ee.d_zeta = ptemp.zeta - ptempin.zeta
                ee.d_delta = ptemp.delta - ptempin.delta

    @classmethod
    def from_sixinput(cls, sixinput, classes=elements):
        other_info = {}

        line_data, rest, iconv = _expand_struct(sixinput, convert=classes)

        ele_names = [dd[0] for dd in line_data]
        elements = [dd[2] for dd in line_data]

        line = cls(elements=elements, element_names=ele_names)

        other_info['rest'] = rest
        other_info['iconv'] = iconv

        return line, other_info

    @classmethod
    def from_madx_sequence(cls,
                           sequence,
                           classes=elements,
                           ignored_madtypes=[],
                           exact_drift=False):

        line = cls(elements=[], element_names=[])

        return _from_madx_sequence(line, sequence, classes, ignored_madtypes, exact_drift)
