import pickle,time,gzip

import numpy as np

import pysixtrack

def savedat(out,fname):
  fmt="%5d"*2 +" %26.15e"*13 +'\n'
  fh=gzip.open(fname,'w')
  for iel,p in enumerate(ptrack):
    for partid in range(npart):
      fh.write(fmt%(iel,partid,
                    p.x[partid],p.px[partid],
                    p.y[partid],p.py[partid],
                    p.tau[partid],p.pt[partid],p.delta[partid],
                    p.chi,p.qratio,p.mratio,
                    p.e0,p.p0c,p.m0))
  fh.close()

# load machine description
sps=pickle.load(open('sps.pickle','rb'))

# particle starting conditions
npart=15
pstart=pysixtrack.Bunch(
         x =np.linspace(0,5e-3,npart),
         px=np.zeros(npart),
         y =np.linspace(0,-3e-3,npart),
         py=np.zeros(npart),
         tau=0.74*np.ones(npart),
         pt=np.zeros(npart),
         e0=26.01692438e9, m0=0.93827205e9)

# element-by-element tracking
ptrack=[]; p=pstart.copy()
for iel,el in enumerate(sps.elems):
    ptrack.append(p.copy())
    el.track(p)

savedat(ptrack,'elem-by-elem.dat.gz')


# many turn tracking
p=pysixtrack.Bunch(x=np.linspace(0,5e-3,100),px=0.,
                   y=np.linspace(0,-3e-3,100),py=0.,
                   tau=0.74,pt=0.,
                   e0=26.01692438e9, m0=0.93827205e9)

ttrack=[]; p=pstart.copy()
before=time.time()
for i in range(25):
  ttrack.append(p.copy())
  sps.track(p)
  now=time.time()
  print(i,"%4.2f sec"%(now-before))
  before=now

savedat(ttrack,'turn-by-turn.dat.gz')
