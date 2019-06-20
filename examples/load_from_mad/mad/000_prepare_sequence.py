import os
from cpymad.madx import Madx

os.system('gfortran headonslice.f -o headonslice')

mad=Madx()
mad.globals.mylhcbeam = 1
mad.globals.on_bb_switch = 1
mad.call('ts_collisions_ats30_en20_IMO380_C7_X160_I1.2_62.3100_60.3200.mask')
mad.input('save, sequence=lhcb1,lhcb2, beam=true, file=lhcwbb.seq;')

