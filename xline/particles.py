import xpart as xp
from .mathlibs import MathlibDefault

class Particles(xp.Particles):

    _m = MathlibDefault

    #def __init__(self, *args, **kwargs):
    #    super().__init__(*args, **kwargs)
    #    self.hide_lost_particles()

    def remove_lost_particles(self):
        self.reorganize()
