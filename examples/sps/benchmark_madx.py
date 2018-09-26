from pyoptics import madlang,optics

#see sps/madx/a001_track_thin.madx
mad=madlang.open('madx/SPS_Q20_thin.seq')
mad.acta_31637.volt=4.5
mad.acta_31637.lag=0.5

import pysixtrack
elems,rest,iconv=mad.sps.expand_struct(pysixtrack.element_types)

pbench=optics.open('madx/track.obs0001.p0001')
sps=pysixtrack.Line(elements= [e[2] for e in elems])

def get_part(pbench,ii):
    pstart=[pbench[n][ii] for n in 'x px y py t pt'.split()]
    pstart=dict(zip('x px y py tau ptau'.split(),pstart))
    prun=pysixtrack.Particles(energy0=pbench.e[ii]*1e9,**pstart)
    return prun

def compare(prun,pbench):
    out=[]
    for att in 'x px y py tau ptau'.split():
        vrun=getattr(prun,att)
        vbench=getattr(pbench,att)
        diff=vrun-vbench
        out.append(abs(diff))
        print(f"{att:<5} {vrun:22.13e} {vbench:22.13e} {diff:22.13g}")
    print(f"max {max(out):21.12e}")
    return max(out)

prun=get_part(pbench,0)
for turn in range(1,30):
    sps.track(prun)
    compare(prun,get_part(pbench,turn))




for ii in range(1,len(iconv)):
    jja=iconv[ii-1]
    jjb=iconv[ii]
    prun=pysixtrack.Particles(


def check_el(p):
    madt=optics.open('madx/track.obs0001.p0001')
    out1=[madt[l][0] for l in 'x px y py y t pt'.split()]
    out2=[getattr(p,l) for l in 'x px y py y tau pt'.split()]
    diff=0
    for a,b in zip(out1,out2):
        diff+=(a-b)**2
        print("%24.17e %24.17e %24.17e"%(a,b,a-b))
    print(np.sqrt(diff))


def trackn(n):
  p=pysixtrack.Particles(x=1e-3,px=0.,y=-0.5e-3,py=0.,tau=0.74,pt=0.,
                     e0=26.01692438e9, m0=0.93827205e9)
  for iel,el in enumerate(sps.elems[:n]):
    print(iel,el)
    el.track(p)
    print("%12.9f %12.9f %12.9f %12.9f %12.9f"%(p.s,p.x,p.px,p.tau,p.pt))
    if abs(p.x)>1:
        break
  print((out[n-1][0]))
  return p



p=trackn(27113)
check_el(p)


def track_turn(n):
  p=pysixtrack.Particles(x=1e-3,px=np.zeros(5000),
                     y=-0.5e-3,py=np.zeros(5000),
                     tau=0.74,pt=np.zeros(5000),
                     e0=26.01692438e9, m0=0.93827205e9)
  out=[]
  before=time.time()
  for i in range(n):
    out.append(p.copy())
    sps.track(p)
    now=time.time()
    print((i,now-before))
    before=now
  return out

tt=track_turn(10)
