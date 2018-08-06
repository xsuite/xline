import sixtracktools
import pysixtrack

six = sixtracktools.SixTrackInput('.')
line, rest, iconv = six.expand_struct(convert=pysixtrack.element_types)
sixdump = sixtracktools.SixDump101('res/dump3.dat')[::2]


print("")
for ii in range(1,5):
    jja=iconv[ii-1]
    jjb=iconv[ii]
    prun=pysixtrack.Particles(**sixdump[ii-1].get_minimal_beam())
    print(f"pysixtr {jja}, x={prun.x}, px={prun.px}")
    for jj in range(jja+1, jjb+1):
        label,elem_type,elem=line[jj]
        elem.track(prun)
        print(f"{label},{elem_type},{str(elem)[:50]}")
        print(f"pysixtr {jj}, x={prun.x}, px={prun.px}")
    pbench=pysixtrack.Particles(**sixdump[ii].get_minimal_beam())
    print(f"sixdump {ii}, x={pbench.x}, px={pbench.px}")
    print("")



