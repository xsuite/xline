import pysixtrack as pyst


for el in pyst.element_list:
    p = pyst.Particles()
    el().track(p)
