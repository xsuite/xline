import numpy as np

from .base_classes import Element
from . import elements
from .particles import Particles

from .loader_sixtrack import _expand_struct
from .loader_mad import _from_madx_sequence


class Line(Element):
    _description = [
        ("elements", "", "List of elements", ()),
        ("element_names", "", "List of element names", ()),
    ]
    _extra = []

    def __len__(self):
        assert len(self.elements) == len(self.element_names)
        return len(self.elements)

    def to_dict(self, keepextra=False):
        out = {}
        out["elements"] = [el.to_dict(keepextra) for el in self.elements]
        out["element_names"] = self.element_names[:]
        return out

    @classmethod
    def from_dict(cls, dct, keepextra=True):
        self = cls(elements=[], element_names=[])
        for el in dct["elements"]:
            eltype = getattr(elements, el["__class__"])
            newel = eltype.from_dict(el)
            self.elements.append(newel)
        self.element_names = dct["element_names"]
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
                newele = getattr(elements, type_name)(**ee._asdict())
                self.elements.append(newele)

        # Append the names
        self.element_names += line.element_names

        assert len(self.elements) == len(self.element_names)

        return self

    def track(self, p):
        ret = None
        for el in self.elements:
            ret = el.track(p)
            if ret is not None:
                break
        return ret

    def track_elem_by_elem(self, p, start=True, end=False):
        out = []
        if start:
            out.append(p.copy())
        for el in self.elements:
            ret = el.track(p)
            if ret is not None:
                break
            out.append(p.copy())
        if end:
            out.append(p.copy())
        return out

    def insert_element(self, idx, element, name):
        self.elements.insert(idx, element)
        self.element_names.insert(idx, name)
        assert len(self.elements) == len(self.element_names)

    def append_element(self, element, name):
        self.elements.append(element)
        self.element_names.append(name)
        assert len(self.elements) == len(self.element_names)

    def get_length(self):
        thick_element_types = (elements.Drift, elements.DriftExact)

        ll = 0
        for ee in self.elements:
            if isinstance(ee, thick_element_types):
                ll += ee.length
        return ll

    def get_s_elements(self, mode="upstream"):
        thick_element_types = (elements.Drift, elements.DriftExact)

        assert mode in ["upstream", "downstream"]
        s_prev = 0
        s = []
        for ee in self.elements:
            if mode == "upstream":
                s.append(s_prev)
            if isinstance(ee, thick_element_types):
                s_prev += ee.length
            if mode == "downstream":
                s.append(s_prev)
        return s

    def remove_inactive_multipoles(self, inplace=False):
        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if isinstance(ee, (elements.Multipole)):
                aux = [ee.hxl, ee.hyl] + ee.knl + ee.ksl
                if np.sum(np.abs(np.array(aux))) == 0.0:
                    continue
            newline.append_element(ee, nn)

        if inplace:
            self.elements.clear()
            self.element_names.clear()
            self.append_line(newline)
            return self
        else:
            return newline

    def remove_zero_length_drifts(self, inplace=False):
        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if isinstance(ee, (elements.Drift, elements.DriftExact)):
                if ee.length == 0.0:
                    continue
            newline.append_element(ee, nn)

        if inplace:
            self.elements.clear()
            self.element_names.clear()
            self.append_line(newline)
            return self
        else:
            return newline

    def merge_consecutive_drifts(self, inplace=False):
        newline = Line(elements=[], element_names=[])

        for ee, nn in zip(self.elements, self.element_names):
            if len(newline.elements) == 0:
                newline.append_element(ee, nn)
                continue

            if isinstance(ee, (elements.Drift, elements.DriftExact)):
                prev_ee = newline.elements[-1]
                prev_nn = newline.element_names[-1]
                if isinstance(prev_ee, (elements.Drift, elements.DriftExact)):
                    prev_ee.length += ee.length
                    prev_nn += nn
                else:
                    newline.append_element(ee, nn)
            else:
                newline.append_element(ee, nn)

        if inplace:
            self.elements.clear()
            self.element_names.clear()
            self.append_line(newline)
            return self
        else:
            return newline

    def get_elements_of_type(self, types):
        if not hasattr(types, "__iter__"):
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

    def find_closed_orbit(
        self, p0c, guess=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], method="Nelder-Mead"
    ):
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
                [pcl.x, pcl.px, pcl.y, pcl.py, pcl.sigma, pcl.delta]
            )

            return coord_out

        def _CO_error(coord):
            return np.sum((_one_turn_map(coord) - coord) ** 2)

        if method == "get_guess":
            res = type("", (), {})()
            res.x = guess
        else:
            import scipy.optimize as so

            res = so.minimize(
                _CO_error, np.array(guess), tol=1e-20, method=method
            )

        pcl = Particles(p0c=p0c)

        pcl.x = res.x[0]
        pcl.px = res.x[1]
        pcl.y = res.x[2]
        pcl.py = res.x[3]
        pcl.zeta = res.x[4]
        pcl.delta = res.x[5]

        return pcl

    def enable_beambeam(self):

        for ee in self.elements:
            if isinstance(ee, (elements.BeamBeam4D, elements.BeamBeam6D)):
                ee.enabled = True

    def disable_beambeam(self):

        for ee in self.elements:
            if isinstance(ee, (elements.BeamBeam4D, elements.BeamBeam6D)):
                ee.enabled = False

    def beambeam_store_closed_orbit_and_dipolar_kicks(
        self,
        particle_on_CO,
        separation_given_wrt_closed_orbit_4D=True,
        separation_given_wrt_closed_orbit_6D=True,
    ):

        self.disable_beambeam()
        closed_orbit = self.track_elem_by_elem(particle_on_CO)

        self.enable_beambeam()

        for ie, ee in enumerate(self.elements):
            # to transfer to beambeam.py

            if ee.__class__.__name__ == "BeamBeam4D":
                if separation_given_wrt_closed_orbit_4D:
                    ee.x_bb += closed_orbit[ie].x
                    ee.y_bb += closed_orbit[ie].y

                # Evaluate dipolar kick
                ptemp = closed_orbit[ie].copy()
                ptempin = ptemp.copy()

                ee.track(ptemp)

                ee.d_px = ptemp.px - ptempin.px
                ee.d_py = ptemp.py - ptempin.py

            elif ee.__class__.__name__ == "BeamBeam6D":
                if not separation_given_wrt_closed_orbit_6D:
                    raise ValueError("Not implemented!")

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

        other_info["rest"] = rest
        other_info["iconv"] = iconv

        return line, other_info

    @classmethod
    def from_madx_sequence(
        cls, sequence, classes=elements, ignored_madtypes=[], exact_drift=False
    ):

        line = cls(elements=[], element_names=[])

        return _from_madx_sequence(
            line, sequence, classes, ignored_madtypes, exact_drift
        )

    # error handling (alignment, multipole orders, ...):

    def add_offset_error_to(self, element, dx=0, dy=0):
        # will raise error if element not present:
        idx_el = self.elements.index(element)
        el_name = self.element_names[idx_el]
        if not dx and not dy:
            return
        xyshift = elements.XYShift(dx=dx, dy=dy)
        inv_xyshift = elements.XYShift(dx=-dx, dy=-dy)
        self.insert_element(idx_el, xyshift, el_name + "_offset_in")
        self.insert_element(idx_el + 2, inv_xyshift, el_name + "_offset_out")

    def add_tilt_error_to(self, element, angle):
        # will raise error if element not present:
        idx_el = self.elements.index(element)
        el_name = self.element_names[idx_el]
        if not angle:
            return
        srot = elements.SRotation(angle=angle)
        inv_srot = elements.SRotation(angle=-angle)
        self.insert_element(idx_el, srot, el_name + "_tilt_in")
        self.insert_element(idx_el + 2, inv_srot, el_name + "_tilt_out")

    def add_multipole_error_to(self, element, knl=[], ksl=[]):
        # will raise error if element not present:
        idx_el = self.elements.index(element)
        el_name = self.element_names[idx_el]
        # normal components
        knl = np.trim_zeros(knl, trim="b")
        if len(element.knl) < len(knl):
            element.knl += [0] * (len(knl) - len(element.knl))
        for i, component in enumerate(knl):
            element.knl[i] += component
        # skew components
        ksl = np.trim_zeros(ksl, trim="b")
        if len(element.ksl) < len(ksl):
            element.ksl += [0] * (len(ksl) - len(element.ksl))
        for i, component in enumerate(ksl):
            element.ksl[i] += component

    def apply_madx_errors(self, error_table):
        """Applies MAD-X error_table (with multipole errors,
        dx and dy offset errors and dpsi tilt errors)
        to existing elements in this Line instance.

        Return error_table names which were not found in the
        elements of this Line instance (and thus not treated).

        Example via cpymad:
            madx = cpymad.madx.Madx()

            # (...set up lattice and errors in cpymad...)

            seq = madx.sequence.some_lattice
            # store already applied errors:
            madx.command.esave(file='lattice_errors.err')
            madx.command.readtable(
                file='lattice_errors.err', table="errors")
            errors = madx.table.errors

            pysixtrack_line, _ = Line.from_madx_sequence(seq)
            pysixtrack_line.apply_madx_errors(errors)
        """
        max_multipole_err = 0
        # check for errors in table which cannot be treated yet:
        for error_type in error_table.keys():
            if error_type == "name":
                continue
            if any(error_table[error_type]):
                if error_type in ["dx", "dy", "dpsi"]:
                    # available alignment error
                    continue
                elif error_type[:1] == "k" and error_type[-1:] == "l":
                    # available multipole error
                    order = int("".join(c for c in error_type if c.isdigit()))
                    max_multipole_err = max(max_multipole_err, order)
                else:
                    print(
                        f'Warning: MAD-X error type "{error_type}"'
                        " not implemented yet."
                    )

        elements_not_found = []
        for i_line, element_name in enumerate(error_table["name"]):
            if not element_name in self.element_names:
                elements_not_found.append(element_name)
                continue
            element = self.elements[self.element_names.index(element_name)]

            # add offset
            try:
                dx = error_table["dx"][i_line]
            except KeyError:
                dx = 0
            try:
                dy = error_table["dy"][i_line]
            except KeyError:
                dy = 0
            self.add_offset_error_to(element, dx, dy)

            # add tilt
            try:
                dpsi = error_table["dpsi"][i_line]
                self.add_tilt_error_to(element, angle=dpsi * 180 / np.pi)
            except KeyError:
                pass

            # add multipole error
            knl = [
                error_table[f"k{o}l"][i_line]
                for o in range(max_multipole_err + 1)
            ]
            ksl = [
                error_table[f"k{o}sl"][i_line]
                for o in range(max_multipole_err + 1)
            ]
            self.add_multipole_error_to(element, knl, ksl)

        return elements_not_found
