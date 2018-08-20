import sixtracktools
import pysixtrack
import pysixtrack.helpers as hp


import numpy as np

# Read sixtrack input
six = sixtracktools.SixTrackInput('.')
p0c_eV = six.initialconditions[-3]*1e6

# Build pysixtrack line
line, rest, iconv = six.expand_struct(convert=pysixtrack.element_types)

# Disable BB elements
ind_BB4D, namelistBB4D, listBB4D = hp.get_elems_of_type(line, 'BeamBeam4D')
for bb in listBB4D: bb.enabled = False

# Find closed orbit
ring = hp.Ring(line, p0c=p0c_eV)
closed_orbit = ring.find_closed_orbit()

# Re-enable beam-beam
for bb in listBB4D: bb.enabled = True

# Add closed orbit to separation (as assumed in sixtrack)
for bb, ibb in zip(listBB4D, ind_BB4D):
    bb.Delta_x+=closed_orbit[ibb].x
    bb.Delta_y+=closed_orbit[ibb].y

# Evaluate kick at CO location
for bb, ibb in zip(listBB4D, ind_BB4D):

    ptemp = closed_orbit[ibb].copy()
    ptempin = ptemp.copy()

    bb.track(ptemp)

    Dpx = ptemp.px - ptempin.px
    Dpy = ptemp.py - ptempin.py

    bb.Dpx_sub = Dpx
    bb.Dpy_sub = Dpy






# Compare closed orbit against sixtrack
# Load sixtrack tracking data
sixdump_all = sixtracktools.SixDump101('res/dump3.dat')
# Assume first particle to be on the closed orbit
Nele_st = len(iconv)
sixdump_CO = sixdump_all[::2][:Nele_st]
x_CO = np.array([pp.x for pp in closed_orbit])
y_CO = np.array([pp.y for pp in closed_orbit])
x_CO_at_st_ele = x_CO[iconv]
y_CO_at_st_ele = y_CO[iconv]
print('Max C.O. discrepancy in x %.2e m'%np.max(np.abs(x_CO_at_st_ele-sixdump_CO.x)))
print('Max C.O. discrepancy in y %.2e m'%np.max(np.abs(y_CO_at_st_ele-sixdump_CO.y)))


# Compare tracking results
sixdump = sixdump_all[1::2] # Particle with deviation from CO

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