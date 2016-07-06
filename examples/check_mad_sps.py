from pyoptics import madlang,optics

#see sps/madx/a001_track_thin.madx
mad=madlang.open('sps/madx/SPS_Q20_thin.seq')
mad.acta_31637.volt=4.5
mad.acta_31637.lag=0.5
out,rest=mad.sps.expand_struct()

import pysixtrack

out,rest=mad.sps.expand_struct(pysixtrack.convert)
elems=zip(*out)[1]
sps=pysixtrack.Block(elems)

import numpy as np
npart=10


def check_el(p):
    madt=optics.open('sps/madx/track.obs0002.p0001')
    out1=[madt[l][0] for l in 'x px y py y t pt e'.split()]
    out1[7]*=1e9
    out2=[getattr(p,l) for l in 'x px y py y tau pt E'.split()]
    diff=0
    for a,b in zip(out1,out2):
        diff+=(a-b)**2
        print "%24.17e %24.17e %24.17e"%(a,b,a-b)
    print np.sqrt(diff)

def trackn(n):
  p=pysixtrack.Bunch(x=1e-3,px=0.,y=-0.5e-3,py=0.,tau=0.74,pt=0.,
                     e0=26.01692438e9, m0=0.93827205e9)
  for iel,el in enumerate(sps.elems[:n]):
    print iel,el
    el.track(p)
    print "%12.9f %12.9f %12.9f %12.9f"%(p.s,p.x,p.px,p.tau)
  print(out[n-1][0])
  return p


p=trackn(27)
check_el(p)



p=Bunch(p0c=26e9,px=0,pt=0)
sps.track(p)


