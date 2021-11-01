import xpart as xp
from .mathlibs import MathlibDefault

class XlineTestParticles(xp.Particles):

    _m = MathlibDefault

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.hide_lost_particles()
        self.delta.flags.writeable=False
        self.rvv.flags.writeable=False
        self.rpp.flags.writeable=False
        self.psigma.flags.writeable=False

    def remove_lost_particles(self):
        self.reorganize()
