import os
from cpymad.madx import Madx

os.system('gfortran headonslice.f -o headonslice')

mad=Madx()
mad.call('ts_collisions_ats30_en20_IMO380_C7_X160_I1.2_62.3100_60.3200.mask')
