from pyoptics import madlang, optics
import pysixtrack

# see sps/madx/a001_track_thin.madx
mad = madlang.open("madx/SPS_Q20_thin.seq")
mad.acta_31637.volt = 4.5
mad.acta_31637.lag = 0.5

elems, rest, iconv = mad.sps.expand_struct(pysixtrack.element_types)

pbench = optics.open("madx/track.obs0001.p0001")
sps = pysixtrack.Line(elements=[e[2] for e in elems])


def get_part(pbench, ii):
    pstart = [pbench[n][ii] for n in "x px y py t pt".split()]
    pstart = dict(zip("x px y py tau ptau".split(), pstart))
    prun = pysixtrack.Particles(energy0=pbench.e[ii] * 1e9, **pstart)
    return prun


def compare(prun, pbench):
    out = []
    for att in "x px y py tau ptau".split():
        vrun = getattr(prun, att)
        vbench = getattr(pbench, att)
        diff = vrun - vbench
        out.append(abs(diff))
        print(f"{att:<5} {vrun:22.13e} {vbench:22.13e} {diff:22.13g}")
    print(f"max {max(out):21.12e}")
    return max(out)


prun = get_part(pbench, 0)
for turn in range(1, 30):
    sps.track(prun)
    compare(prun, get_part(pbench, turn))
