import numpy as np

from .base_classes import Element
from .beambeam import BeamBeam4D
from .beambeam import BeamBeam6D

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


class Drift(Element):
    """Drift in expanded form"""
    _description = [
        ("length", "m", "Length of the drift", 0),
    ]

    def track(self, p):
        length = self.length
        rpp = p.rpp
        xp = p.px * rpp
        yp = p.py * rpp
        p.x += xp * length
        p.y += yp * length
        p.zeta += length * (p.rvv - (1 + (xp**2 + yp**2) / 2))
        p.s += length


class DriftExact(Drift):
    """Drift in exact form"""
    _description = [
        ("length", "m", "Length of the drift", 0),
    ]

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
    """ Multipole """
    _description = [
        ('knl', 'm^-n', "Normalized integrated strength of normal components", 0),
        ('ksl', 'm^-n', "Normalized integrated strength of skew components", 0),
        ('hxl', 'rad', "Rotation angle of the reference trajectory in the horizzontal plane", 0),
        ('hyl', 'rad', "Rotation angle of the reference trajectory in the vertical plane", 0),
        ('length', 'm', "Length of the orginating thick multipole", 0),
    ]

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


class Cavity(Element):
    """Radio-frequency cavity"""
    _description = [
        ('voltage', 'V', "Integrated energy change", 0),
        ('frequency', 'Hz', "Frequency of the cavity", 0),
        ('lag', 'degree', "Delay in the cavity sin(lag - w tau)", 0),
    ]

    def track(self, p):
        sin = p._m.sin
        pi = p._m.pi
        k = 2 * pi * self.frequency / p.clight
        tau = p.zeta / p.rvv / p.beta0
        phase = self.lag * pi / 180 - k * tau
        p.add_to_energy(p.chi * self.voltage * sin(phase))



class XYShift(Element):
    """shift of the reference"""
    _description = [
        ('dx', 'm', 'Horizontal shift', 0),
        ('dy', 'm', 'Vertical shift', 0)
    ]

    def track(self, p):
        p.x -= self.dx
        p.y -= self.dy


class SRotation(Element):
    """anti-clockwise rotation of the reference frame"""
    _description = [
        ('angle', '', 'Rotation angle', 0),
    ]

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



