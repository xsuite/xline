import time
import gzip
import math

from pysixtrack import (
    Drift,
    DriftExact,
    Multipole,
    Align,
    Cavity,
    Block,
    Particles,
)


def savedata(out, fname):
    fmt = "%5d" * 2 + " %26.15e" * 13 + "\n"
    fh = gzip.open(fname, "w")
    for iel, p in enumerate(ptrack):
        fh.write(
            fmt
            % (
                iel,
                0,
                p.x,
                p.px,
                p.y,
                p.py,
                p.tau,
                p.pt,
                p.delta,
                p.chi,
                p.qratio,
                p.mratio,
                p.e0,
                p.p0c,
                p.m0,
            )
        )
    fh.close()


# Load machine description
sps = eval(open("sps.dat").read())


# particle starting conditions

pstart = Particles(
    x=5e-3,
    px=0,
    y=-3e-3,
    py=0,
    tau=0.74,
    pt=0,
    e0=26.01692438e9,
    m0=0.93827205e9,
    mathlib=math,
)

# element-by-element tracking

ptrack = []
p = pstart.copy()
for iel, el in enumerate(sps.elems):
    ptrack.append(p.copy())
    el.track(p)

savedata(ptrack, "elem-by-elem2.dat.gz")


# many turn tracking

ttrack = []
p = pstart.copy()
for i in range(5):
    ttrack.append(p.copy())
    before = time.time()
    for j in range(100):
        sps.track(p)
    now = time.time()
    print((i + 1) * 100, "%4.2f sec" % (now - before))

savedata(ttrack, "turn-by-turn2.dat.gz")
