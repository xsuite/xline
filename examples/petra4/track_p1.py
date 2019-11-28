from cpymad.madx import Madx

import pysixtrack


class Madout:
    def __init__(self):
        self.out = []

    def write(self, ll):
        self.out.append(ll)

    def __repr__(self):
        return "".join(ll.decode() for ll in self.out)


madout = Madout()
mad = Madx(stdout=madout)

mad.call("h7ba_n8.seq")
mad.beam(energy=6, sequence="ring", particle="electron", radiate=False)
mad.use("ring")
mad.twiss()
print(mad.table.summ.q1, mad.table.summ.q2)

nslices = 4
mad.select(flag="makethin", clear=True)
mad.select(flag="makethin", class_="sbend", slice=nslices)
mad.select(flag="makethin", class_="quadrupole", slice=nslices)
mad.makethin(sequence="ring")
mad.use(sequence="ring")
print(mad.table.summ.q1, mad.table.summ.q2)
twiss = mad.twiss()
print(mad.table.summ.q1, mad.table.summ.q2)


twissout = pysixtrack.Particles.from_madx_twiss(
    mad.twiss(betx=1, bety=1, x=0.001)
)

line = pysixtrack.Line.from_madx_sequence(mad.sequence.ring)
part = pysixtrack.Particles()
part.x = 0.001
pysixout = pysixtrack.Particles.from_list(
    line.track_elem_by_elem(part, start=False, end=True)
)


def mkd(name, t1, t2, ii, jj):
    v1 = getattr(t1, name)[ii]
    v2 = getattr(t2, name)[jj]
    print(f"{name:4}   {v2:20} {v2:20} {v2-v1:20}")
    return v2 - v1


for ii in range(len(twissout.s)):
    sm = twissout.s[ii]
    sp = pysixout.s[ii]
    print(
        f"{ii:3} {twiss.name[ii]:20}"
        "{line.element_names[ii]:20} {line.elements[ii]}"
    )
    res = mkd("s", twissout, pysixout, ii, ii)
    if abs(res) > 1e-6:
        break
