# Dofulator

This library provides a degrees of freedom (DoF) calculator for molecular dynamics
simulations which use rigid geometry constraints.
Per-atom and per-direction DoF are calculated using the methods described
in [J. Chem. Theory Comput. 2024, 20, 23, 10615–10624](https://doi.org/10.1021/acs.jctc.4c00957),
which are applicable to both rigid bodies and semi-rigid fragments.

The core library is provided as a C API (see [the usage docs](doc/usage.md#C) or [dofulator.h](src/dofulator.h)).
A Python wrapper is also provided, and includes some plugins compatible with
[MDAnalysis](https://github.com/MDAnalysis/mdanalysis) (some example usage is shown [below](#mdanalysis-examples)).

For further details, see [the documentation](doc/README.md).


## Installation

See [here](doc/installation.md) for full installation details.

To quickly set up for usage with MDAnalysis:
```
pip install "dofulator[mdanalysis]@git+https://github.com/CTCMS-UQ/dofulator.git"
```

NOTE: This will require a compatible BLAS and LAPACK implementation to be installed.
MKL and OpenBLAS have both been confirmed to work.

### Linux

Linux users can install BLAS and LAPACK via their package manager, or load the relevant module
if using an HPC environment.
Note, depending on your Linux distribution you may need to install the cblas and lapacke versions, possibly with a -dev or -devel suffix.

### Mac

Mac users can install OpenBLAS via homebrew.
The following flags (or similar) will likely be needed with the above command to properly locate the installation:
```
--config-settings=cmake.define.BLA_VENDOR=OpenBLAS --config-settings=cmake.define.CMAKE_PREFIX_PATH=/opt/homebrew/Cellar/openblas/<VERSION_NUMBER>/ --config-settings=cmake.define.CMAKE_C_FLAGS='-I/opt/homebrew/Cellar/openblas/<VERSION_NUMBER/include'
```
Make sure to set the correct version number for the version you have installed.

### Windows

For Windows users, BLAS and LAPACK can be installed via Conda with:
```
TODO
```
Additionally, note that the MSVC compiler is *NOT* currently supported, since it does not
fully comply with modern C standards.
Instead, the build defaults to the Zig C compiler bundled in the `ziglang` python package.

## MDAnalysis examples

### Directional DoF of all atoms

```python
import MDAnalysis as mda
from dofulator import MDADofulator
u = mda.Universe('topology.tpr', 'trajectory.trr')
d = MDADofulator(
    u.atoms,            # Include all relevant atoms in the dofulator context
    rigid_bonds=[b for b in u.bonds if b.type == '3'],          # Set bond type 3 as rigid
    rigid_bodies=[r for r in u.residues if r.resname == 'H2O'], # Water treated as rigid bodies
    mode='directional'  # Options are 'atomic' (default) or 'directional'
)
d.run()
# Access results via d.results
```

### Temperature profile across system

For local temperature calculation, the `LocalTemperature` class makes
use of MDAnalysis' `AtomGroup`s, which can dynamically update as the trajectory plays.
```python
import MDAnalysis as mda
import numpy as np
from dofulator import MDADofulator, LocalTemperature
u = mda.Universe('topology.tpr', 'trajectory.trr')
zrange = np.linspace(0, u.dimensions[2], 11) # 10 bins in the z direction
t = LocalTemperature(
    # List of atom selections, each of which will get temperature calculated on each frame
    [u.select_atoms(f'prop z >= {zrange[i]} and prop z < {zrange[i+1]}', updating=True)
        for i in range(len(zrange) - 1)],

    # Dofulator context for calculating DoF in each selection
    MDADofulator(
        u.atoms,            # Include all atoms which could be selected
        rigid_bonds=[b for b in u.bonds if b[0].type == 'H' or b[1].type == 'H'], # Rigid bonds to hydrogens
    )
)
t.run()
# t.results contains an n_frames x 10 array with temperature of each bin on each frame
```

## Citation

When using dofulator in published work, please cite the following paper:
1. S. Sanderson, S. R. Tee, and D. J. Searles; Local Temperature Measurement in Molecular Dynamics Simulations with Rigid Constraints *Journal of Chemical Theory and Computation* **2024** *20* (23), 10615–10624. DOI: [10.1021/acs.jctc.4c00957](https://doi.org/10.1021/acs.jctc.4c00957)
