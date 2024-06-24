
from .c_dofulator import CDofulator as Dofulator

try:
    # Will fail if MDAnalysis not available
    from .mda_dofulator import MDADofulator
    from .local_temperature import LocalTemperature
except ImportError:
    pass

