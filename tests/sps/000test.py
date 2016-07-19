import pickle,time

import pysixtrack

import numpy as np

# load machine description
sps=pickle.load(open('sps.pickle','rb'))


# element-by-element tracking
p=pysixtrack.Bunch(x=np.linspace(0,5e-3,100),px=0.,
                   y=np.linspace(0,-3e-3,100),py=0.,
                   tau=0.74,pt=0.,
                   e0=26.01692438e9, m0=0.93827205e9)

ptrack=[]
for iel,el in enumerate(sps.elems):
    ptrack.append(p.copy())
    el.track(p)

# single turn tracking
p=pysixtrack.Bunch(x=np.linspace(0,5e-3,100),px=0.,
                   y=np.linspace(0,-3e-3,100),py=0.,
                   tau=0.74,pt=0.,
                   e0=26.01692438e9, m0=0.93827205e9)
sps.track(p)
pend=p.copy()


# many turn tracking
p=pysixtrack.Bunch(x=np.linspace(0,5e-3,100),px=0.,
                   y=np.linspace(0,-3e-3,100),py=0.,
                   tau=0.74,pt=0.,
                   e0=26.01692438e9, m0=0.93827205e9)

out=[]
before=time.time()
for i in range(25):
  out.append(p.copy())
  sps.track(p)
  now=time.time()
  print(i,"%4.2f sec"%(now-before))
  before=now

