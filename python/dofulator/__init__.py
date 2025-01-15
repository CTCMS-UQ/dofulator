try:
    # Add scipy-openblas32 library directory to path in case it was used to
    # build so it can be loaded by c_dofulator
    import scipy_openblas32
    import sys
    sys.path.append(scipy_openblas32.get_lib_dir())
except ImportError:
    pass

from .c_dofulator import CDofulator as Dofulator

try:
    # Will fail if MDAnalysis not available
    from .mda_dofulator import MDADofulator
    from .local_temperature import LocalTemperature
except ImportError:
    pass

