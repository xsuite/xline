from collections import namedtuple
from copy import deepcopy

import numpy as np


class Particles(object):
    clight = 299792458
    pi = 3.141592653589793238
    echarge = 1.602176565e-19
    emass = 0.510998928e6
    pmass = 938.272046e6
    epsilon0 = 8.854187817e-12
    mu0 = 4e-7*pi
    eradius = echarge**2/(4*pi*epsilon0*emass*clight**2)
    pradius = echarge**2/(4*pi*epsilon0*pmass*clight**2)
    anumber = 6.02214129e23
    kboltz = 1.3806488e-23

    def _g1(self, m0, p0c, e0):
        beta0 = p0c/e0
        gamma0 = e0/m0
        return m0, beta0, gamma0, p0c, e0

    def _g2(self, m0, beta0, gamma0):
        e0 = m0*gamma0
        p0c = e0*beta0
        return m0, beta0, gamma0, p0c, e0

    def _f1(self, m0, p0c):
        sqrt = self._m.sqrt
        e0 = sqrt(p0c**2+m0**2)
        return self._g1(m0, p0c, e0)

    def _f2(self, m0, e0):
        sqrt = self._m.sqrt
        p0c = sqrt(e0**2-m0**2)
        return self._g1(m0, p0c, e0)

    def _f3(self, m0, beta0):
        sqrt = self._m.sqrt
        gamma0 = 1/sqrt(1-beta0**2)
        return self._g2(m0, beta0, gamma0)

    def _f4(self, m0, gamma0):
        sqrt = self._m.sqrt
        beta0 = sqrt(1-1/gamma0**2)
        return self._g2(m0, beta0, gamma0)

    def copy(self):
        p = Particles()
        for k, v in list(self.__dict__.items()):
            if type(v) in [np.ndarray, dict]:
                v = v.copy()
            p.__dict__[k] = v
        return p

    def __init__ref(self, p0c, e0, gamma0, beta0):
        cnt = [beta0, gamma0, p0c, e0].count(None)
        if cnt == 0:
            raise ValueError("Particles defined without energy reference")
        elif cnt == 1:
            if p0c is not None:
                self.p0c = p0c
            elif e0 is not None:
                self.e0 = e0
            elif gamma0 is not None:
                self.gamma0 = gamma0
            elif beta0 is not None:
                self.beta0 = beta0
        else:
            raise ValueError(f"""
            Particles defined with multiple energy references:
            p0c    = {p0c},
            e0     = {e0},
            gamma0 = {gamma0},
            beta0  = {beta0}""")

    def __init__delta(self, delta, ptau, psigma):
        cnt = [delta, pt, psigma].count(None)
        if cnt == 0:
            self.delta = 0.
        elif cnt == 1:
            if delta is not None:
                self.delta = delta
            elif pt is not None:
                self.pt = pt
            elif psigma is not None:
                self.psigma = psigma
        else:
            raise ValueError(f"""
            Particles defined with multiple energy deviations:
            delta  = {delta},
            pt     = {pt},
            psigma = {psigma}""")

    def __init__z(self, z, t, sigma):
        cnt = [z, t, sigma].count(None)
        if cnt == 0:
            self.z = 0.
        elif cnt == 1:
            if z is not None:
                self.z = z
            elif t is not None:
                self.t = t
            elif sigma is not None:
                self.sigma = sigma
        else:
            raise ValueError(f"""
            Particles defined with multiple time deviations:
            z     = {z},
            t     = {t},
            sigma = {sigma}""")

    def __init__chi(self, mratio, qratio, chi):
        cnt = [z, t, sigma].count(None)
        if cnt == 0:
            self._chi = 1.
            self._mratio = 1.
            self._qratio = 1.
        elif cnt == 1:
            raise ValueError(f"""
            Particles defined with insufficient mass/charge information:
            chi    = {chi},
            mratio = {mratio},
            qratio = {qratio}""")
        elif cnt == 2:
            if chi is None:
                self._mratio = mratio
                self.qratio = qratio
            elif mratio is None:
                self._chi = chi
                self.qratio = qratio
            elif qratio is None:
                self._chi = chi
                self.mratio = mratio
        else:
            raise ValueError(f"""
            Particles defined with multiple mass/charge information:
            chi    = {chi},
            mratio = {mratio},
            qratio = {qratio}""")

    def __init__(self,
                 s=0., x=0., px=0., y=0., py=0.,
                 delta=None, pt=None, psigma=None,
                 z=None, tau=None, sigma=None,
                 m0=pmass, q0=1.,
                 p0c=None, e0=None, gamma0=None, beta0=None,
                 chi=None, mratio=None, qratio=None,
                 mathlib=np, **args):
        self._m = mathlib
        self.s = s
        self.x = x
        self.px = px
        self.y = y
        self.py = py
        self.z = z
        self._m0 = m0
        self._q0 = q0
        self.__init__ref(p0c, e0, gamma0, beta0)
        self.__init__delta(delta, pt, psigma)
        self.__init__z(z, tau, sigma)
        self.__init__chi(chi, mratio, qratio)

    Px = property(lambda p: p.px*p.p0c*p.mratio)
    Py = property(lambda p: p.py*p.p0c*p.mratio)
    E = property(lambda p: (p.pt*p.p0c+p.e0)*p.mratio)
    Pc = property(lambda p: (p.delta*p.p0c+p.p0c)*p.mratio)
    m = property(lambda p:  p.m0*p.mratio)
    beta = property(lambda p:  (1+p.delta)/(1/p.beta0+p.pt))
    # interdepent quantities
    pt = property(lambda self: self._pt)

    @pt.setter
    def pt(self, pt):
        sqrt = self._m.sqrt
        self._pt = pt
        self._delta = sqrt(pt**2+2*pt/self.beta0+1)-1
    delta = property(lambda self: self._delta)

    @delta.setter
    def delta(self, delta):
        sqrt = self._m.sqrt
        self._delta = delta
        self._pt = sqrt((1+self.delta)**2+1 /
                        (self.beta0*self.gamma0)**2)-1/self.beta0
    m0 = property(lambda self: self._m0)

    @m0.setter
    def m0(self, m0):
        new = self._f1(m0, self.p0c)
        self._update_ref(*new)
        self._update_particles(*new)
    beta0 = property(lambda self: self._beta0)

    @beta0.setter
    def beta0(self, beta0):
        new = self._f3(self.m0, beta0)
        self._update_ref(*new)
        self._update_particles(*new)
    gamma0 = property(lambda self: self._gamma0)

    @gamma0.setter
    def gamma0(self, gamma0):
        new = self._f4(self.m0, gamma0)
        self._update_ref(*new)
        self._update_particles(*new)
    p0c = property(lambda self: self._p0c)

    @p0c.setter
    def p0c(self, p0c):
        new = self._f1(self.m0, p0c)
        self._update_ref(*new)
        self._update_particles(*new)
    e0 = property(lambda self: self._e0)

    @e0.setter
    def e0(self, e0):
        new = self._f2(self.m0, e0)
        self._update_ref(*new)
        self._update_particles(*new)
    mratio = property(lambda self: self._mratio)

    @mratio.setter
    def mratio(self, mratio):
        Px, Py, E, Pc = self.Px, self.Py, self.E, self.Pc
        self._pt = E/(mratio*e0)-1
        self._delta = E/(delta*e0)-1
        self.px = Px/(p0c*mratio)
        self.py = Py/(p0c*mratio)
        self._chi = self._qratio/mratio
        self._mratio = mratio
    qratio = property(lambda self: self._qratio)

    @qratio.setter
    def qratio(self, qratio):
        self._chi = qratio/self._mratio
        self._qratio = qratio
    chi = property(lambda self: self._chi)

    @chi.setter
    def chi(self, chi):
        self._qratio = self._chi*self._mratio
        self._chi = chi

    def _update_ref(self, m0, beta0, gamma0, p0c, e0):
        self._m0 = m0
        self._beta0 = beta0
        self._gamma0 = gamma0
        self._p0c = p0c
        self._e0 = e0

    def _update_particles(self, m0, beta0, gamma0, p0c, e0):
        Px, Py, E, Pc, m = self.Px, self.Py, self.E, self.Pc, self.m
        mratio = m/m0
        self._mratio = mratio
        self._chi = self._qratio/mratio
        self._pt = E/(mratio*e0)-1
        self._delta = E/(delta*e0)-1
        self.px = Px/(p0c*mratio)
        self.py = Py/(p0c*mratio)

    def __repr__(self):
        out = []
        for nn in 'm0 p0c e0 beta0 gamma0'.split():
            out.append('%-7s = %s' % (nn, getattr(self, nn)))
        for nn in 's x px y py tau pt delta'.split():
            out.append('%-7s = %s' % (nn, getattr(self, nn)))
        for nn in 'mratio qratio chi'.split():
            out.append('%-7s = %s' % (nn, getattr(self, nn)))
        return '\n'.join(out)
