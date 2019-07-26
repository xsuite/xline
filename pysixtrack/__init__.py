__version__ = "0.0.1"

from . import base_classes
from . import elements

from .line import Line

from .particles import Particles

from .be_beamfields import BB6Ddata
from .loader_mad import MadPoint

elements.Line = Line

element_list = [
    elements.Drift,
    elements.DriftExact,
    elements.Multipole,
    elements.Cavity,
    elements.XYShift,
    elements.SRotation,
    elements.RFMultipole,
    elements.BeamMonitor,
    elements.DipoleEdge,
    elements.Line,
]
