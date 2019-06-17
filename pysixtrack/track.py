import numpy as np
import math

from copy import deepcopy

from scipy.constants import e as qe
from scipy.constants import c as clight
from .be_beambeam.gaussian_fields import get_Ex_Ey_Gx_Gy_gauss
from .particles import Particles
from .be_beambeam import BB6D
from .be_beambeam import BB6Ddata

_factorial = np.array([1,
                       1,
                       2,
                       6,
                       24,
                       120,
                       720,
                       5040,
                       40320,
                       362880,
                       3628800,
                       39916800,
                       479001600,
                       6227020800,
                       87178291200,
                       1307674368000,
                       20922789888000,
                       355687428096000,
                       6402373705728000,
                       121645100408832000,
                       2432902008176640000])


class Element(object):
    def __init__(self, **nargs):
        not_allowed = set(nargs.keys()) - set(self.__slots__)
        if len(not_allowed) > 0:
            cname = self.__class__.__name__
            raise NameError(f"{not_allowed} not allowed by {cname}")
        for name, default in zip(self.__slots__, self.__defaults__):
            setattr(self, name, deepcopy(nargs.get(name, default)))

    def _slots(self):
        names = self.__slots__
        return zip(names, [getattr(self, name) for name in names])

    def __repr__(self):
        args = [f"{name}={val}" for name, val in self._slots()]
        args = ", ".join(args)
        return f"{self.__class__.__name__}({args})"

    def _asdict(self):
        return dict(self._slots())


class Drift(Element):
    __slots__ = ('length',)
    __units__ = ('meter',)
    __defaults__ = (0,)

    def track(self, p):
        length = self.length
        rpp = p.rpp
        xp = p.px * rpp
        yp = p.py * rpp
        p.x += xp * length
        p.y += yp * length
        p.zeta += length * (p.rvv - (1 + (xp**2 + yp**2) / 2))
        p.s += length


class DriftExact(Element):
    __slots__ = ('length',)
    __units__ = ('meter',)
    __defaults__ = (0,)

    def track(self, p):
        sqrt = p._m.sqrt
        length = self.length
        opd = 1 + p.delta
        lpzi = length / sqrt(opd**2 - p.px**2 - p.py**2)
        p.x += p.px * lpzi
        p.y += p.py * lpzi
        p.zeta += p.rvv * length - opd * lpzi
        p.s += length


class Multipole(Element):
    __slots__ = ('knl', 'ksl', 'hxl', 'hyl', 'length')
    __units__ = ('meter^-order', 'meter^-order', 'radian', 'radian', 'meter')
    __defaults__ = ([], [], 0, 0, 0)

    @property
    def order(self):
        return max(len(self.knl), len(self.ksl)) - 1

    def track(self, p):
        order = self.order
        length = self.length
        knl = np.array(self.knl)
        ksl = np.array(self.ksl)
        if len(knl) < len(ksl):
            nknl = np.zeros(order + 1, dtype=knl.dtype)
            nknl[:len(knl)] = knl
            knl = nknl
        elif len(knl) > len(ksl):
            nksl = np.zeros(order + 1, dtype=ksl.dtype)
            nksl[:len(ksl)] = ksl
            ksl = nksl
        x = p.x
        y = p.y
        chi = p.chi
        dpx = knl[order]
        dpy = ksl[order]
        for ii in range(order, 0, -1):
            zre = (dpx * x - dpy * y) / ii
            zim = (dpx * y + dpy * x) / ii
            dpx = knl[ii - 1] + zre
            dpy = ksl[ii - 1] + zim
        dpx = -chi * dpx
        dpy = chi * dpy
        # curvature effect kick
        hxl = self.hxl
        hyl = self.hyl
        delta = p.delta
        if (hxl != 0 or hyl != 0):
            b1l = chi * knl[0]
            a1l = chi * ksl[0]
            hxlx = hxl * x
            hyly = hyl * y
            if (length > 0):
                hxx = hxlx / length
                hyy = hyly / length
            else:  # non physical weak focusing disabled (SixTrack mode)
                hxx = 0
                hyy = 0
            dpx += hxl + hxl * delta - b1l * hxx
            dpy -= hyl + hyl * delta - a1l * hyy
            p.zeta -= chi * (hxlx - hyly)
        p.px += dpx
        p.py += dpy


class XYShift(Element):
    """shift of the reference"""
    __slots__ = ('dx', 'dy')
    __units__ = ('meter', 'meter')
    __defaults__ = (0, 0)

    def track(self, p):
        p.x -= self.dx
        p.y -= self.dy


class SRotation(Element):
    """anti-clockwise rotation of the reference frame"""
    __slots__ = ('angle',)
    __units__ = ('degree',)
    __defaults__ = (0,)

    def track(self, p):
        deg2rag = p._m.pi / 180
        cz = p._m.cos(self.angle * deg2rag)
        sz = p._m.sin(self.angle * deg2rag)
        xn = cz * p.x + sz * p.y
        yn = -sz * p.x + cz * p.y
        p.x = xn
        p.y = yn
        xn = cz * p.px + sz * p.py
        yn = -sz * p.px + cz * p.py
        p.px = xn
        p.py = yn


class Cavity(Element):
    __slots__ = ('voltage', 'frequency', 'lag')
    __units__ = ('volt', 'hertz', 'degree')
    __defaults__ = (0, 0, 0)

    def track(self, p):
        sin = p._m.sin
        pi = p._m.pi
        k = 2 * pi * self.frequency / p.clight
        tau = p.zeta / p.rvv / p.beta0
        phase = self.lag * pi / 180 - k * tau
        p.add_to_energy(p.chi * self.voltage * sin(phase))


class RFMultipole(Element):
    __slots__ = ('voltage', 'frequency', 'knl', 'ksl', 'pn', 'ps')
    __units__ = ('volt', 'hertz', [], [], [], [])
    __defaults__ = (0, 0, 0, 0, 0)


class BeamMonitor(Element):
    __slots__ = ('num_stores', 'start', 'skip', 'max_particle_id',
                 'min_particle_id', 'is_rolling', 'is_turn_ordered')
    __units__ = ('', 'turn', 'turn', '', '', '', '')
    __defaults__ = (0, 0, 1, 0, 0, False, True)

    def offset(self, particle):
        _offset = -1
        nn = self.max_particle_id >= self.min_particle_id \
            and (self.max_particle_id - self.min_particle_id + 1) or -1
        assert(self.is_turn_ordered)

        if particle.turn >= self.start and nn > 0 and \
                particle.partid >= self.min_particle_id and \
                particle.partid <= self.max_particle_id:
            turns_since_start = particle.turns - self.start
            store_index = turns_since_start // self.skip
            if store_index < self.num_stores:
                pass
            elif self.is_rolling:
                store_index = store_index % self.num_stores
            else:
                store_index = -1

            if store_index >= 0:
                _offset = store_index * nn + particle.partid

        return _offset

    def track(self, particle):
        pass


class Line(Element):
    __slots__ = ('elements', 'element_names' )
    __defaults__ = ([], [])

    def track(self, p):
        for el in self.elements:
            el.track(p)

    def track_elem_by_elem(self, p):
        out = []
        for el in self.elements:
            out.append(p.copy())
            el.track(p)
        return out

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

    @classmethod
    def fromline(cls, line):
        return cls().append_line(line)

    def find_closed_orbit(self, p0c, guess=[0.,0.,0.,0.,0.,0.],
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

        CO = self.track_elem_by_elem(pcl)

        return CO


class Monitor(Element):
    __slots__ = ('data',)
    __defaults__ = ([],)

    def track(self, p):
        self.data.append(p.copy)


class BeamBeam4D(Element):
    __slots__ = (
        'charge',
        'sigma_x',
        'sigma_y',
        'beta_r',
        'min_sigma_diff',
        'x_bb',
        'y_bb',
        'd_px',
        'd_py',
        'enabled')
    __units__ = ('e', 'm', 'm', [],
                 'm', 'm', 'm', [], [], [])
    __defaults__ = (0., 0., 0., 0.,
                    1e-28, 0., 0., 0., 0., True)

    def track(self, p):
        if self.enabled:
            charge = p.qratio * p.q0
            x = p.x - self.x_bb
            px = p.px
            y = p.y - self.y_bb
            py = p.py

            chi = p.chi

            beta = p.beta0 / p.rvv
            p0c = p.p0c * qe

            Ex, Ey = get_Ex_Ey_Gx_Gy_gauss(
                x, y, self.sigma_x, self.sigma_y, min_sigma_diff=1e-10,
                skip_Gs=True, mathlib=p._m)

            fact_kick = chi * self.charge * qe * charge * qe\
                * (1. + beta * self.beta_r) / (p0c * (beta + self.beta_r))
            
            px += (fact_kick * Ex - self.d_px)
            py += (fact_kick * Ey - self.d_py)

            p.px = px
            p.py = py


class BeamBeam6D(Element):
    __slots__ = ([
        'phi',
        'alpha',
        'x_bb_co',
        'y_bb_co',
        'charge_slices',
        'zeta_slices',
        'sigma_11',
        'sigma_12',
        'sigma_13',
        'sigma_14',
        'sigma_22',
        'sigma_23',
        'sigma_24',
        'sigma_33',
        'sigma_34',
        'sigma_44',
        'x_co',
        'px_co',
        'y_co',
        'py_co',
        'zeta_co',
        'delta_co',
        'd_x',
        'd_px',
        'd_y',
        'd_py',
        'd_zeta',
        'd_delta',
        'min_sigma_diff',
        'threshold_singular',
        'enabled'])

    __units__ = tuple(len(__slots__) * [[]])
    __defaults__ = tuple((len(__slots__)-3) * [0.] + [1e-28, 1e-28, True])

    def track(self, p):
        if self.enabled:
            bb6data = BB6Ddata.BB6D_init(
                qe,
                self.phi,
                self.alpha,
                self.x_bb_co,
                self.y_bb_co,
                self.charge_slices,
                self.zeta_slices,
                self.sigma_11,
                self.sigma_12,
                self.sigma_13,
                self.sigma_14,
                self.sigma_22,
                self.sigma_23,
                self.sigma_24,
                self.sigma_33,
                self.sigma_34,
                self.sigma_44,
                self.x_co,
                self.px_co,
                self.y_co,
                self.py_co,
                self.zeta_co,
                self.delta_co,
                self.min_sigma_diff,
                self.threshold_singular,
                self.d_x,
                self.d_px,
                self.d_y,
                self.d_py,
                self.d_zeta,
                self.d_delta,
                self.enabled)
            x_ret, px_ret, y_ret, py_ret, zeta_ret, delta_ret = BB6D.BB6D_track(
                p.x, p.px, p.y, p.py, p.zeta, p.delta, p.q0 * qe,
                p.p0c / clight * qe, bb6data)
            self._last_bb6data = bb6data
            p.x = x_ret
            p.px = px_ret
            p.y = y_ret
            p.py = py_ret
            p.zeta = zeta_ret
            p.delta = delta_ret


classes = [cls for cls in globals().values() if isinstance(cls, type)]
elements = [cls for cls in classes if issubclass(cls, Element)]
__all__ = [cls.__name__ for cls in elements]
__all__.append('element_types')

element_types = dict((cls.__name__, cls) for cls in elements)
