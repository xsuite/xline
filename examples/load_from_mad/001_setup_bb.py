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
   
    bb_elements = []
    bb_names = []
    bb_xyz_points = []
    for ee in seq.elements:
        if ee.base_type.name == 'beambeam':
            eename = ee.name
            bb_elements.append(ee)
            bb_names.append(eename)
            bb_xyz_points.append(MadPoint(eename+':1', mad))
    
    return bb_elements, bb_names, bb_xyz_points

def add_beambeam(mad,sequence,name,s,from_):
    el=mad.command.beambeam.clone(name)
    mad.seqedit(sequence=sequence)
    mad.install(element=name,at=s,from_=from_)
    mad.endedit()
    return mad.elements[name]

from cpymad.madx import Madx
import pysixtrack

mad=Madx()
mad.options.echo=False;
mad.options.warn=False;
mad.options.info=False;

mad.call('mad/lhcwbb.seq')


bb_ele_b1, bb_names_b1, bb_xyz_b1 = get_bb_names_and_xyz_points(mad, seq_name='lhcb1')

# Install BB in B2
for bb1 in bb_ele_b1:
    mad.seqedit(sequence='lhcb2')
    newbb=mad.command.beambeam.clone(bb1.name,
            sigx=1e-3, sigy=1e-3, xma=0., yma=0, charge=0., 
            comments=bb1.comments)
    mad.install(element=bb1.name, at=bb1.at,from_=bb1['from'])
    mad.endedit()

#bb_ele_b1, bb_names_b1, bb_xyz_b1 = get_bb_names_and_xyz_points(mad, seq_name='lhcb1')
bb_ele_b2, bb_names_b2, bb_xyz_b2 = get_bb_names_and_xyz_points(mad, seq_name='lhcb2')

# A check
assert len(bb_names_b1)==len(bb_names_b2)
for nbb1, nbb2 in zip(bb_names_b1, bb_names_b2):
    assert(nbb1==nbb2)


i_bb = 10
pb1 = bb_xyz_b1[i_bb]
pb2 = bb_xyz_b2[i_bb]







#line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.lhcb1)
