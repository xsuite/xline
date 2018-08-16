import sixtracktools
import pysixtrack

import numpy as np

six = sixtracktools.SixTrackInput('.')
line, rest, iconv = six.expand_struct(convert=pysixtrack.element_types)
sixdump_all = sixtracktools.SixDump101('res/dump3.dat')

N_st_ele = len(iconv) #Number of elements in SixTrack
N_pyst_ele = len(line) #Number of elements in PySixTrack

# extract closed orbit and test particle
sixdump_CO = sixdump_all[::2]
sixdump = sixdump_all[1::2]

# Type of PyST elements
pyst_ele_type_list = [ele[1] for ele in line]
pyst_used_types = set(pyst_ele_type_list)

# Index of bb elements in PyST
pyst_ind_BB4D = np.where([ss=='BeamBeam4D' for ss in pyst_ele_type_list])[0]

# Build list of bb elements in PyST
bb_ele_list = [line[ind] for ind in pyst_ind_BB4D]

# Index of bb elem in ST
st_ind_BB4D = np.array([iconv.index(ind) for ind in pyst_ind_BB4D])

#Find closed-orbit at beam-beam interactions
st_CO_exit_BB4D = sixdump_CO[st_ind_BB4D]

# Adjust delta definitions
for ibb, bb_ele in enumerate(bb_ele_list):
    bb = bb_ele[2]

    delta_x_st = bb.Delta_x
    bb.Delta_x = st_CO_exit_BB4D.x[ibb]+delta_x_st

    delta_y_st = bb.Delta_y
    bb.Delta_y = st_CO_exit_BB4D.y[ibb]+delta_y_st

# Evaluate kick at CO location
for ibb, bb_ele in enumerate(bb_ele_list):
    bb = bb_ele[2]
    ptemp = pysixtrack.Particles(**st_CO_exit_BB4D[ibb].get_minimal_beam())
    ptempin = ptemp.copy()

    bb.track(ptemp)

    Dpx = ptemp.px - ptempin.px
    Dpy = ptemp.py - ptempin.py



def compare(prun,pbench):
    out=[]
    for att in 'x px y py delta sigma'.split():
        vrun=getattr(prun,att)
        vbench=getattr(pbench,att)
        diff=vrun-vbench
        out.append(abs(diff))
        print(f"{att:<5} {vrun:22.13e} {vbench:22.13e} {diff:22.13g}")
    print(f"max {max(out):21.12e}")
    return max(out)

print("")
for ii in range(1,len(iconv)):
    jja=iconv[ii-1]
    jjb=iconv[ii]
    prun=pysixtrack.Particles(**sixdump[ii-1].get_minimal_beam())
    print(f"\n-----sixtrack={ii} sixtracklib={jja} --------------")
    #print(f"pysixtr {jja}, x={prun.x}, px={prun.px}")
    for jj in range(jja+1, jjb+1):
        label,elem_type,elem=line[jj]
        elem.track(prun)
        print(f"{jj} {label},{str(elem)[:50]}")
    pbench=pysixtrack.Particles(**sixdump[ii].get_minimal_beam())
    #print(f"sixdump {ii}, x={pbench.x}, px={pbench.px}")
    print("-----------------------")
    out=compare(prun,pbench)
    print("-----------------------\n\n")
    if out>1e-13:
        print("Too large discrepancy")
        break


