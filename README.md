# Dofulator

This library provides a degrees of freedom (DoF) calculator for molecular dynamics
simulations which use rigid geometry constraints.
Per-atom and per-direction DoF are calculated using the methods described
in [J. Chem. Theory Comput. 2024, 20, 23, 10615–10624](https://doi.org/10.1021/acs.jctc.4c00957),
which are applicable to both rigid bodies and semi-rigid fragments.

The core library is provided as a C API (see [the docs](doc/usage.md#C) or [dofulator.h](src/dofulator.h)).
For further details, see [the documentation](doc/README.md).

A Python wrapper is also provided, and includes some plugins compatible with
[MDAnalysis](https://github.com/MDAnalysis/mdanalysis).
Some example usage is shown [below](#Examples).


## Installation

See [here](doc/installation.md) for full installation details.

To quickly set up for usage with MDAnalysis:
```
pip install git+https://github.com/CTCMS-UQ/dofulator[mdanalysis]
```


## Examples

### Directional DoF of all atoms

```python
import MDAnalysis as mda
from dofulator import MDADofulator
u = mda.Universe('trajectory.traj')
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
u = mda.Universe('trajectory.traj')
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
