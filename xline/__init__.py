__version__ = "0.0.11"

raise ImportError("xline is not anymore supported by xsuite\n"
    "    Instead of `xline.Line` use `xtrack.Line`\n"
    "    Instead of `xline.Particles` use `xpart.Particles`\n"
    "    Instead of `xline.Drift` use `xtrack.Drift`\n"
    "    Instead of `xline.Multipole` use `xtrack.Multipole`\n"
    "    ... etc ...")

from . import base_classes
from . import elements
from .elements import *

from .line import Line

from .particles import XlineTestParticles

from .be_beamfields import BB6Ddata
from .loader_mad import MadPoint, mad_benchmark

