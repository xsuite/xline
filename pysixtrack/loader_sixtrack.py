from math import factorial

import numpy as np

from .base_classes import Element
from . import elements 

clight = 299792458
pi = np.pi

def bn_mad(bn_mad, n, sign):
    return sign*bn_mad*factorial(n-1)


def bn_rel(bn16, bn3, r0, d0, sign):
    out = []
    for nn, (a, b) in enumerate(zip(bn16, bn3)):
        n = nn+1
        sixval = d0*a*b*r0**(1-n)*10**(3*n-6)
        out.append(bn_mad(sixval, n, sign))
    return out

def _expand_struct(self, convert=elements):
    elems = []
    count = {}
    icount = 0
    iconv = []
    names = []
    rest = []
    Drift = convert.Drift
    Multipole = convert.Multipole
    Cavity = convert.Cavity
    XYShift = convert.XYShift
    SRotation = convert.SRotation
    #Line = convert.Line
    BeamBeam4D = convert.BeamBeam4D
    BeamBeam6D = convert.BeamBeam6D
    exclude = False
    # add special elenents
    if 'CAV' in self.iter_struct():
        self.single['CAV'] = [12*self.ition, self.u0, self.harm, 0]
    for nnn in self.iter_struct():
        exclude = False
        ccc = count.setdefault(nnn, 0)
        if len(self.single[nnn]) == 7:
            etype, d1, d2, d3, d4, d5, d6 = self.single[nnn]
        else:
            etype, d1, d2, d3, = self.single[nnn]
            d4, d5, d6 = None, None, None
        elem = None
        if nnn in self.align:
            dx, dy, tilt = self.align[nnn][ccc]
            tilt = tilt*180e-3/pi
            dx *= 1e-3
            dy *= 1e-3
            hasshift = abs(dx)+abs(dy) > 0
            hastilt = abs(tilt) > 0
            if hasshift:
                names.append(nnn+'_preshift')
                elems.append(XYShift(dx=dx, dy=dy))
                icount += 1
            if hastilt:
                names.append(nnn+'_pretilt')
                elems.append(SRotation(angle=tilt))
                icount += 1
        if etype in [0, 25]:
            elem = Drift(length=d3)
            if d3 > 0:
                exclude = True
        elif abs(etype) in [1, 2, 3, 4, 5, 7, 8, 9, 10]:
            bn_six = d1
            nn = abs(etype)
            sign = -etype/nn
            madval = bn_mad(bn_six, nn, sign)
            knl = [0]*(nn-1)+[madval]
            ksl = [0]*nn
            if sign == 1:
                knl, ksl = ksl, knl
            elem = Multipole(knl=knl, ksl=ksl, hxl=0, hyl=0, length=0)
        elif etype == 11:
            knl, ksl = self.get_knl(nnn, ccc)
            hxl = 0
            hyl = 0
            l = 0
            # beaware of the case of thick bend
            # see beambeam example where mbw has the length
            if d3 == -1:
                hxl = -d1
                l = d2
                knl[0] = hxl
            elif d3 == -2:
                hyl = -d1  # strange sign!!!
                l = d2
                ksl[0] = hyl
            elem = Multipole(knl=knl, ksl=ksl, hxl=hxl, hyl=hyl, length=l)
        elif etype == 12:
            # e0=self.initialconditions[-1]
            # p0c=np.sqrt(e0**2-self.pma**2)
            # beta0=p0c/e0
            v = d1*1e6
            freq = d2*clight/self.tlen
            # print(v,freq)
            elem = Cavity(voltage=v, frequency=freq, lag=180-d3)
        elif etype == 20:
            thisbb = self.bbelements[nnn]
            if type(thisbb) is self.classes['BeamBeam4D']:
                elem = BeamBeam4D(**thisbb._asdict())
            elif type(thisbb) is self.classes['BeamBeam6D']:
                elem = BeamBeam6D(**thisbb._asdict())
            else:
                raise ValueError('What?!')
        else:
            rest.append([nnn]+self.single[nnn])
        if elem is not None:
            elems.append(elem)
            names.append(nnn)
        if nnn in self.align:
            if hastilt:
                names.append(nnn+'_posttilt')
                elems.append(SRotation(angle=-tilt))
                icount += 1
            if hasshift:
                names.append(nnn+'_postshift')
                elems.append(XYShift(dx=-dx, dy=-dy))
                icount += 1
        if elem is not None:
            if not exclude:
                iconv.append(icount)
            icount += 1
        count[nnn] = ccc+1
    #newelems = [dict(i._asdict()) for i in elems]
    types = [i.__class__.__name__ for i in elems]
    return list(zip(names, types, elems)), rest, iconv
