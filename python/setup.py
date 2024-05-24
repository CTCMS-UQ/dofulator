from setuptools import setup, Extension
from Cython.Build import cythonize
import os
import sys

# Copied from MDAnalysis to handle numpy build dependency
def get_numpy_include():
    # Obtain the numpy include directory. This logic works across numpy
    # versions.
    # setuptools forgets to unset numpy's setup flag and we get a crippled
    # version of it unless we do it ourselves.
    import builtins

    builtins.__NUMPY_SETUP__ = False
    try:
        import numpy as np
    except ImportError:
        print('*** package "numpy" not found ***')
        print('Dofulator requires a version of NumPy, even for setup.')
        print('Please get it from http://numpy.scipy.org/ or install it through '
              'your package manager.')
        sys.exit(-1)
    return np.get_include()

__module_dir__ = os.path.abspath(os.path.dirname(__file__))
setup(
    include_dirs=[get_numpy_include()],
    ext_modules=cythonize([
        Extension(
            'dofulator.dofulator',
            [os.path.join(__module_dir__, 'dofulator', 'dofulator.pyx')],
            libraries=['dofulator'],
            define_macros=[('NPY_NO_DEPRECATED_API', 'NPY_1_7_API_VERSION')],
        )
    ]),
)

