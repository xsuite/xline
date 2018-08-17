import numpy as np
import math
import numba

from scipy.constants import e as qe

from .gaussian_fields import get_Ex_Ey_Gx_Gy_gauss

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
        not_allowed=set(nargs.keys())-set(self.__slots__)
        if len(not_allowed)>0:
            cname=self.__class__.__name__
            raise NameError(f"{not_allowed} not allowed by {cname}")
        for name, default in zip(self.__slots__, self.__defaults__):
            setattr(self, name, nargs.get(name, default))

    def _slots(self):
        names = self.__slots__
        return zip(names, [getattr(self, name) for name in names])

    def __repr__(self):
        args = [f"{name}={val}" for name, val in self._slots()]
        args = ", ".join(args)
        return f"{self.__class__.__name__}({args})"

    def as_dict(self):
        return dict(self._slots())


class Drift(Element):
    __slots__ = ('length',)
    __units__ = ('meter',)
    __defaults__ = (0,)

    def track(self, p):
        length = self.length
        rpp = p.rpp
        xp = p.px*rpp
        yp = p.py*rpp
        p.x += xp*length
        p.y += yp*length
        p.zeta += length*(p.rvv-(1+(xp**2+yp**2)/2))
        p.s += length


class DriftExact(Element):
    __slots__ = ('length',)
    __units__ = ('meter',)
    __defaults__ = (0,)

    def track(self, p):
        sqrt = p._m.sqrt
        length = self.length
        rpp = p.rpp
        xp = p.px*rpp
        yp = p.py*rpp
        p.x += xp*length
        p.y += yp*length
        p.zeta += length*(p.rvv+(xp**2+yp**2)/2)
        p.s += length


class Multipole(Element):
    __slots__ = ('knl', 'ksl', 'hxl', 'hyl', 'length')
    __units__ = ('meter^-order', 'meter^-order', 'radian', 'radian', 'meter')
    __defaults__ = ([], [], 0, 0, 0)

    @property
    def order(self):
        return max(len(self.knl), len(self.ksl))-1

    def track(self, p):
        order = self.order
        length = self.length
        knl = np.array(self.knl)
        ksl = np.array(self.ksl)
        if len(knl) < len(ksl):
            nknl = np.zeros(order+1, dtype=knl.dtype)
            nknl[:len(knl)] = knl
            knl = nknl
        elif len(knl) > len(ksl):
            nksl = np.zeros(order+1, dtype=ksl.dtype)
            nksl[:len(ksl)] = ksl
            ksl = nksl
        x = p.x
        y = p.y
        chi = p.chi
        dpx = knl[order]
        dpy = ksl[order]
        for ii in range(order, 0, -1):
            zre = (dpx*x-dpy*y)/ii
            zim = (dpx*y+dpy*x)/ii
            dpx = knl[ii-1]+zre
            dpy = ksl[ii-1]+zim
        dpx = -chi*dpx
        dpy = chi*dpy
        # curvature effect kick
        hxl = self.hxl
        hyl = self.hyl
        delta = p.delta
        if (length > 0):
            b1l = chi*knl[0]
            a1l = chi*ksl[0]
            hxx = hxl/length*x
            hyy = hyl/length*y
            dpx += hxl + hxl*delta - b1l*hxx
            dpy -= hyl + hyl*delta - a1l*hyy
            p.zeta -= chi*(hxx-hyy)*length
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
        deg2rag = p._m.pi/180
        cz = p._m.cos(self.angle*deg2rag)
        sz = p._m.sin(self.angle*deg2rag)
        print(cz, sz)
        xn = cz*p.x + sz*p.y
        yn = -sz*p.x + cz*p.y
        p.x = xn
        p.y = yn
        xn = cz*p.px + sz*p.py
        yn = -sz*p.px + cz*p.py
        p.px = xn
        p.py = yn


class Cavity(Element):
    __slots__ = ('voltage', 'frequency', 'lag')
    __units__ = ('volt', 'hertz', 'degree')
    __defaults__ = (0, 0, 0)

    def track(self, p):
        sin = p._m.sin
        pi = p._m.pi
        k = 2*pi*self.frequency/p.clight
        tau = p.zeta/p.rvv/p.beta0
        phase = self.lag*pi/180-k*tau
        p.add_to_energy(p.chi*self.voltage*sin(phase))


class RFMultipole(Element):
    __slots__ = ('voltage', 'frequency', 'knl', 'ksl', 'pn', 'ps')
    __units__ = ('volt', 'hertz', [], [], [], [])
    __defaults__ = (0, 0, 0, 0, 0)


class Line(Element):
    __slots__ = ('elements',)
    __defaults__ = ([],)

    def track(self, p):
        for el in self.elements:
            el.track(p)

    def track_elem_by_elem(self,p):
        out=[]
        for el in self.elements:
            out.append(p.copy())
            el.track(p)
        return out

class BeamBeam4D(Element):
    __slots__ = ('q_part', 'N_part', 'sigma_x', 'sigma_y', 'beta_s',
                 'min_sigma_diff', 'Delta_x', 'Delta_y', 'Dpx_sub', 'Dpy_sub')
    __units__ = tuple(len(__slots__)*[[]])
    __defaults__ = tuple(len(__slots__)*[0.])

    def track(self, p):
        charge = p.qratio*p.q0*qe
        x = p.x - self.Delta_x
        px = p.px
        y = p.y - self.Delta_y
        py = p.py

        chi = p.chi

        beta = p.beta0/p.rvv
        p0c = p.p0c*qe

        Ex, Ey = get_Ex_Ey_Gx_Gy_gauss(x, y, self.sigma_x, self.sigma_y,
                                       min_sigma_diff=1e-10, skip_Gs=True, mathlib=p._m)

        fact_kick = chi * self.N_part * self.q_part * charge * \
            (1. + beta * self.beta_s)/(p0c*(beta + self.beta_s))

        px += (fact_kick*Ex - self.Dpx_sub)
        py += (fact_kick*Ey - self.Dpy_sub)

        p.px = px
        p.py = py


class BeamBeam6D(Element):
    __slots__ = tuple(('ibsix xang xplane h_sep v_sep ' +
                       'sigma_xx sigma_xxp sigma_xpxp sigma_yy sigma_yyp ' +
                       'sigma_ypyp sigma_xy sigma_xyp sigma_xpy sigma_xpyp strengthratio').split())
    __units__ = tuple(len(__slots__)*[[]])
    __defaults__ = tuple(len(__slots__)*[0.])


classes = [cls for cls in globals().values() if isinstance(cls, type)]
elements = [cls for cls in classes if issubclass(cls, Element)]
__all__ = [cls.__name__ for cls in elements]
__all__.append('element_types')

element_types = dict((cls.__name__, cls) for cls in elements)
