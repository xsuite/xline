import numpy as np
from numpy.random import uniform
import matplotlib.pyplot as plt

import pysixtrack

p = pysixtrack.Particles()

rect_aperture = pysixtrack.elements.LimitRect(
    min_x=-1e-2, max_x=2e-2, min_y=-0.5e-2, max_y=2.5e-2
)

ellip_aperture = pysixtrack.elements.LimitEllipse(a=1.5e-2, b=0.3e-2)

# Test scalar
assert rect_aperture.track(p) is None
assert ellip_aperture.track(p) is None

p.x = 100.0
assert rect_aperture.track(p) == 0
assert ellip_aperture.track(p) == 0

# Test vector

N_part = 10000

p.x = uniform(low=-3e-2, high=3e-2, size=N_part)
p.y = uniform(low=-3e-2, high=3e-2, size=N_part)

p.state = np.ones_like(p.x, dtype=np.int)

rect_aperture.track(p)
ellip_aperture.track(p)


plt.close("all")
fig1 = plt.figure(1)
ax = fig1.add_subplot(111)

ax.plot(p.x, p.y, ".b")
ax.plot(p.lost_particles[0].x, p.lost_particles[0].y, "r.")
ax.plot(p.lost_particles[1].x, p.lost_particles[1].y, "g.")

plt.show()
