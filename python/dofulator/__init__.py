
from .dofulator import CDofulator as Dofulator

try:
    # Will fail if MDAnalysis not available
    from .mda_dofulator import MDADofulator
except ImportError:
    pass
