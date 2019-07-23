import pickle

import numpy as np
import matplotlib.pylab as plt

from cpymad.madx import Madx

import pysixtrack
from pysixtrack.particles import Particles

from sc_controller import SC_controller

p0c = 25.92e9
intensity=2e11 
neps_x=2e-6
neps_y=2e-6
dpp_rms=1.5e-3
bunchlength_rms = 0.22
V_RF_MV = 4.5
lag__RFdeg = 180.
max_distance = 35. #6.9 #25.

mad = Madx()
mad.options.echo=False
mad.options.info=False
mad.warn=False 

mad.call('madx/sps_thin.seq')
mad.use('sps')
def install_spacecharge(mad, seq_name, name, s, mode='Bunched'): 

    mad.input("""
        ??NAME?? : beambeam, sigx   = 1., sigy    = 1.,
            xma    = 0.,
            yma    = 0.,
            charge = 0.;
        ??NAME??, slot_id = ??SID??;
        seqedit, sequence=??SEQNAME??;
        install,element=??NAME??,at=??POS??;
        flatten;
        endedit;
        use, sequence=??SEQNAME??;
        """.replace('??NAME??', name).replace(
            '??POS??', '%.10e'%s).replace(
                '??SEQNAME??',seq_name).replace(
                    '??SID??', {'Coasting':'1', 'Bunched':'2'}[mode]))

install_spacecharge(mad, 'sps', 'sc0', s=20.01, mode='Coasting')
install_spacecharge(mad, 'sps', 'sc1', s=20.02)


line, other = pysixtrack.Line.from_madx_sequence(mad.sequence.sps)                             
			     
		             
