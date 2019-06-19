import numpy as np

class MadPoint(object):
    def __init__(self,name,mad):
        self.name=name
        twiss=mad.table.twiss
        survey=mad.table.survey
        idx=np.where(survey.name==name)[0][0]
        self.tx=twiss.x[idx]
        self.ty=twiss.y[idx]
        self.sx=survey.x[idx]
        self.sy=survey.y[idx]
        self.sz=survey.z[idx]
        theta=survey.theta[idx]
        phi=survey.phi[idx]
        psi=survey.psi[idx]
        thetam=np.array([[np.cos(theta) ,           0,np.sin(theta)],
             [          0,           1,         0],
             [-np.sin(theta),           0,np.cos(theta)]])
        phim=np.array([[          1,          0,          0],
            [          0,np.cos(phi)   ,   np.sin(phi)],
            [          0,-np.sin(phi)  ,   np.cos(phi)]])
        psim=np.array([[   np.cos(psi),  -np.sin(psi),          0],
            [   np.sin(psi),   np.cos(psi),          0],
            [          0,          0,          1]])
        wm=np.dot(thetam,np.dot(phim,psim))
        self.ex=np.dot(wm,np.array([1,0,0]))
        self.ey=np.dot(wm,np.array([0,1,0]))
        self.ez=np.dot(wm,np.array([0,0,1]))
        self.sp=np.array([self.sx,self.sy,self.sz])
        self.p=self.sp+ self.ex * self.tx + self.ey * self.ty
    def dist(self,other):
        return np.sqrt(np.sum((self.p-other.p)**2))
    def distxy(self,other):
        dd=self.p-other.p
        return np.dot(dd,self.ex),np.dot(dd,self.ey)

def get_bb_names_and_xyz_points(mad, seq_name):
    mad.use(sequence=seq_name);
    mad.twiss()
    mad.survey()
    
    seq = mad.sequence[seq_name]
   
    bb_names = []
    bb_xyz_points = []
    for ee in seq.elements:
        if ee.base_type.name == 'beambeam':
            eename = ee.name
            bb_names.append(eename)
            bb_xyz_points.append(MadPoint(eename+':1', mad))
    
    return bb_names, bb_xyz_points

def norm(v):
    return np.sqrt(np.sum(v**2))


from cpymad.madx import Madx
import pysixtrack

mad=Madx()
mad.options.echo=False;
mad.options.warn=False;
mad.options.info=False;


# Load sequence
mad.call('mad/lhcwbb.seq')

# Disable beam-beam kicks 
mad.globals.on_bb_charge = 0.

ip_names = [1, 2, 5, 8]

# Retrieve geometry information
mad.use('lhcb1'); mad.twiss(); mad.survey()
IP_xyz_b1 = {}
for ip in ip_names:
    IP_xyz_b1[ip] = MadPoint('ip%d'%ip+':1', mad)

mad.use('lhcb2'); mad.twiss(); mad.survey()
IP_xyz_b2 = {}
for ip in ip_names:
    IP_xyz_b2[ip] = MadPoint('ip%d'%ip+':1', mad)

bb_names_b1, bb_xyz_b1 = get_bb_names_and_xyz_points(
        mad, seq_name='lhcb1')
bb_names_b2, bb_xyz_b2 = get_bb_names_and_xyz_points(
        mad, seq_name='lhcb2')


# A check
assert len(bb_names_b1)==len(bb_names_b2)
for nbb1, nbb2 in zip(bb_names_b1, bb_names_b2):
    assert(nbb1==nbb2.replace('b2_','b1_'))

# # Check that the two reference system are parallel
# assert(norm(pb1.ex-pb2.ex)<1e-10) #1e-4 is a reasonable limit
# assert(norm(pb1.ey-pb2.ey)<1e-10) #1e-4 is a reasonable limit
# assert(norm(pb1.ez-pb2.ez)<1e-10) #1e-4 is a reasonable limit
# 
# # Check that there is no horizontal separation
# assert(np.abs(np.dot(v12, pb1.ez))<1e-10)


import matplotlib.pyplot as plt
plt.close('all')
fig1 = plt.figure(1)

plt.plot(           [pb.p[0] for pb in bb_xyz_b1],
                    [pb.p[2] for pb in bb_xyz_b1], 'b.')
plt.quiver(np.array([pb.p[0] for pb in bb_xyz_b1]),
           np.array([pb.p[2] for pb in bb_xyz_b1]), 
          np.array([pb.ex[0] for pb in bb_xyz_b1]),
          np.array([pb.ex[2] for pb in bb_xyz_b1]))
plt.quiver(np.array([pb.p[0] for pb in bb_xyz_b1]),
           np.array([pb.p[2] for pb in bb_xyz_b1]), 
          np.array([pb.ez[0] for pb in bb_xyz_b1]),
          np.array([pb.ez[2] for pb in bb_xyz_b1]))

plt.plot(           [pb.p[0] for pb in bb_xyz_b2],
                    [pb.p[2] for pb in bb_xyz_b2], 'r.')
plt.quiver(np.array([pb.p[0] for pb in bb_xyz_b2]),
           np.array([pb.p[2] for pb in bb_xyz_b2]), 
          np.array([pb.ex[0] for pb in bb_xyz_b2]),
          np.array([pb.ex[2] for pb in bb_xyz_b2]))
plt.quiver(np.array([pb.p[0] for pb in bb_xyz_b2]),
           np.array([pb.p[2] for pb in bb_xyz_b2]), 
          np.array([pb.ez[0] for pb in bb_xyz_b2]),
          np.array([pb.ez[2] for pb in bb_xyz_b2]))
plt.show()


#line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.lhcb1)
