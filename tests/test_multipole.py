import pysixtrack

p1=pysixtrack.Particles()

p1.x=1 ;p1.y=1
p2=p1.copy()

el1=pysixtrack.elements.RFMultipole(knl=[.5,2,.2],ksl=[.5,3,.1])
el2=pysixtrack.elements.Multipole(knl=el1.knl,ksl=el1.ksl)

el1.track(p1)
el2.track(p2)

assert p1.compare(p2,abs_tol=1e-15)

