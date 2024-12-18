# Using dofulator

* [CLI](#CLI)
* [C](#C)
    - [Setting up the context](#Setting-up-the-context)
    - [Calculating and querying DoF](#Calculating-and-querying-DoF)
    - [Considerations for kinematic loops](#Considerations-for-kinematic-loops)
    - [Querying fragments](#Querying-fragments)
    - [Errors](#Errors)
* [Python](#Python)
* [MDAnalysis](#MDAnalysis)

## CLI

The command-line interface (`dof`) takes an `.xyz` file (optionally with added
bond information) and outputs the total and directional DoF of each atom.
If bond information is not included, it is assumed that all atoms in the file
form a single rigid body.

See the [tests](tests/frontend) directory for example input (`.xyz`) and output
(`.dof`) files, or run `dof --help` for more information.


## C

The library interface works around a dofulator *context*.
Rigid and semi-rigid fragments are defined by adding connections to the context.
Fragments are then finalised, and rigid fragment DoF is pre-calculated from reference geometry.
From this point, repeated calculations can be made with different configurations of the fragments,
and the results queried to get the DoF of atoms of interest.

### Setting up the context
The context is created by calling `dofulator_create(N)`, where `N` is
the total number of atoms in the system:
```C
Dofulator ctx = dofulator_create_context(N);
```
`ctx` is an opaque handle, and should be cleaned up at the end with `dofulator_destroy(ctx)`.


Rigid bodies can be specified by connecting atoms with `dofulator_build_rigid_fragment(ctx, (Bond){i, j})`.
For example, if atoms 0, 1 and 2 form a rigid body, it can be specified with:
```C
dofulator_build_rigid_fragment(ctx, (Bond){0, 1});
dofulator_build_rigid_fragment(ctx, (Bond){0, 2});
```
or alternatively:
```C
dofulator_build_rigid_fragment(ctx, (Bond){0, 1});
dofulator_build_rigid_fragment(ctx, (Bond){1, 2});
```
Note, it doesn't matter which two atoms are provided for an individual call,
only that all atoms in the rigid body are somehow connected to each other by
separate calls.


Semi-rigid fragments are constructed in a similar manner by specifying rigid bonds
between atoms using `dofulator_add_rigid_bond(ctx, (Bond){i, j})`.
For example, if atoms 3, 4 and 5 form a water molecule, with 3 as O and 4 and 5 as H atoms,
then to treat it as having rigid bond lengths but a flexible angle, this would be specified
using:
```C
dofulator_add_rigid_bond(ctx, (Bond){3, 4});
dofulator_add_rigid_bond(ctx, (Bond){3, 5});
```
**IMPORTANT:** Currently, it is not supported for rigid bodies and semi-rigid fragments to
be connected, and attempting to do so will result in an error.
This may be relaxed in future for greater computational efficiency, but for now
rigid sections of a semi-rigid fragment can be achieved by constraining all DoF
in that section.
For example, to treat the water molecule above as rigid, an additional
rigid bond constraining the angle could be added with:
```C
dofulator_add_rigid_bond(ctx, (Bond){4, 5});
```

Once all fragments have been specified, they should be finalised by calling:
```C
dofulator_finalise_fragments(ctx);
```
This prepares the context for DoF calculation, and fragments should not be changed
after this point (i.e. `dofulator_build_rigid_fragment` and `dofulator_add_rigid_bond`
should not be called).
**IMPORTANT:** All functions described below assume that `dofulator_finalise_fragments`
has been called.


If the system has periodic boundaries, these should be specified before any calculations
using `dofulator_set_pbc(ctx, pbc)`, where `pbc` is a tagged union type which should be
treated as follows:
```C
PBC pbc = (PBC){.typ = PBC_NONE}; // Default - no periodicity

// Orthogonal box
pbc = (PBC){.typ = PBC_ORTHO, .lx = x_length, .ly = y_length, .lz = z_length};

// Triclinic box
pbc = (PBC){
    .typ = PBC_TRI,
    .ax = x_length,
    .bx = xy_skew,
    .by = y_length,
    .cx = xz_skew,
    .cy = yz_skew,
    .cz = z_length,
};
```


### Calculating and querying DoF

Since the total atomic DoF of rigid bodies is invariant to translation and rotation,
and the directional components can be determined along an arbitrary axis once the
required matrix has been computed, these values are pre-calculated for
efficiency by calling:
```C
dofulator_precalculate_rigid(ctx, mass, x);
```
where `mass` is a length `N` array of `double`s containing the mass of each
atom (indexed by atom index), and `x` is a flattened `N` by `3` array of atom
positions (i.e. `{x0, y0, z0, x1, y1, z1, ..., xN, yN, zN}`).
The positions should be consistent with the constraint geometry, and are used
to calculate DoF and to define a basis relative to which directional DoF will
be rotated as needed.


For a given system configuration, the DoF of each atom can be determined by
calling
```C
dofulator_calculate(ctx, mass, x);
```
and then queried for a given atom with
```C
double total_atom_dof = dofulator_get_dof_atom(ctx, atom_index);

double directional_dof[3];
dofulator_get_dof_atom_directional(ctx, atom_index, directional_dof);
```
**IMPORTANT:** If `dofulator_build_rigid_fragment` has been called, then
`dofulator_precalculate_rigid` MUST be called before the first call to
`dofulator_calculate`.


If DoF of a list of atoms is required, this can be queried through
```C
AtomTag atoms[] = {0, 2, 4, 3};     // List of atom indices to query
double* dof_totals = dofulator_get_dof_atoms(ctx, 4, atoms, NULL);
// dof_totals = {dof_0, dof_2, dof_4, dof_3}

double* dof_directions = dofulator_get_dof_atoms_directional(ctx, 4, atoms, NULL);
// dof_directions = { \
//   dof_0_x, dof_0_y, dof_0_z, \
//   dof_2_x, dof_2_y, dof_2_z, \
//   dof_4_x, dof_4_y, dof_4_z, \
//   dof_3_x, dof_3_y, dof_3_z \
// }
```
The second argument, `n_atoms`, is the number of atoms in the list.
The third argument can be set to `NULL`, in which case DoF of the first `n_atoms` atoms is returned.
A buffer can be passed in as the final argument in place of `NULL`, in which case the result will
be written to that buffer, and the same pointer will be returned.
Note, it is assumed that the buffer has enough space for `n_atoms` `double`s
for `dofulator_get_dof_atoms`, or for `3 * n_atoms` `double`s for
`dofulator_get_dof_atoms_directional`.
If `NULL` is passed as the buffer (as in the above example), one will be allocated, or `NULL` will be
returned if the allocation fails.


### Considerations for kinematic loops

For topologies with closed kinematic loops, loop closing constraints are handled
by finding the null space of the matrix $`\mathbf{K}`$ which describes motion along those constraints.
This is achieved with singular value decomposition by finding the singular values below
a small fraction of the largest singular value.
With the default null space threshold value of 0.0, the threshold fraction is
determined as `DBL_EPSILON` multiplied by the largest dimension of the $`\mathbf{K}`$
matrix.
This should be accurate for standard cases, but can be changed if needed using
`dofulator_set_null_space_thresh(ctx, thresh)`.
Internally, the absolute value of `thresh` is taken, and is then clamped to a
maximum of 1, so higher values will be interpreted as 1, and negative values
will be interpreted as positive.
The current value can be queried with `thresh = dofulator_get_null_space_thresh(ctx)`.
For some topologies that result in a $`\mathbf{K}`$ matrix where numerical precision
becomes a problem, it may be necessary to raise the cut-off in order to get the
correct total number of DoF, although such cases are very unlikely in
practice.


### Querying fragments

The list of rigid fragments can be iterated over as follows:
```C
FragmentListIter rigid_frags = dofulator_get_rigid_fragments(ctx);
const Fragment* frag;
while ((frag = fragmentlist_iter_next(rigid_frags))) {
    AtomListView atoms = fragment_get_atoms(frag);
    for (size_t i = 0; i < atoms.n_atoms; ++i) {
        AtomTag atom_i = atoms.atoms[i];
        // Do something with the atom index...
    }
}

// OR

for (size_t frag_id = 0; frag_id < rigid_frags.n_fragments; ++frag_id) {
    const Fragment* frag = &rigid_frags.fragments[frag_id];
    // ...
}
```
Note, the former case consumes the iterator, so a new one must be obtained
to repeat the loop.

Similarly, for semi-rigid fragments, a `FragmentListIter` can be obtained
by calling `dofulator_get_semirigid_fragments(ctx)`.


### Errors

Each call to `dofulator_build_rigid_fragment` or `dofulator_add_rigid_bond`
could fail, and hence the return value (of `enum` type `DofulatorResult`)
should be checked. The same is true of `dofulator_finalise_fragments`,
`dofulator_precalculate_rigid`, and `dofulator_calculate`.

Possible values are:

| Value                       | Description |
|-----------------------------|-------------|
| `DOF_SUCCESS`               | No error.   |
| `DOF_UNINITIALISED`         | Received an uninitialised context. |
| `DOF_ALLOC_FAILURE`         | Failed to allocate memory. |
| `DOF_INDEX_ERROR`           | Index out of range. |
| `DOF_MIXED_RIGID_SEMIRIGID` | Tried to link a rigid body to a semi-rigid fragment. |
| `DOF_LAPACK_ERROR`          | Error encountered during a LAPACK call. |




## Python

Dofulator provides both a standalone Python interface, and a wrapper around that which is
compatible with MDAnalysis.
In most cases, the latter (described in the [next section](#MDAnalysis)) will be the simplest to use.
The former is simply a thin layer around the underlying C interface with some convenience functions.

To create a context
```python
from dofulator import Dofulator
import numpy as np

ctx = Dofulator(num_atoms)
```

Rigid bonds which form semi-rigid fragments can be added with
```python
ctx.add_rigid_bond(atom_index_i, atom_index_j)
```

Similarly, atoms can be joined into a rigid body with
```python
ctx.build_rigid_fragment(atom_index_i, atom_index_j)
```

As for the C interface, once all rigid constraints have been added to the context, fragments are finalised
with
```python
ctx.finalise_fragments()
```

Periodic boundary conditions can be disabled with `ctx.set_pbc_none()`, set as orthorhombic with
` ctx.set_pbc_ortho(np.array([lx, ly, lz], dtype=np.float64))`, or set as triclinic with a lower triangular
box matrix from vectors `a`, `b` and `c` using
`ctx.set_pbc_triclinic(np.array([[ax, 0, 0], [bx, by, 0], [cx, cy, cz]], dtype=np.float64))`.

With periodicity set, if any rigid fragments are present (i.e. `ctx.build_rigid_fragment(...)` was called),
then a call must be made to
```python
ctx.precalculate_rigid(masses, positions)
```
where `masses` is a 1D Numpy array of 64 bit floats with length >= the number of atoms (`ctx.n_atoms`),
and `positions` is a 2D Numpy array of 64 bit floats with shape `(ctx.n_atoms, 3)`.

At this point, `ctx.calculate(masses, positions)` can be called repeatedly for various configurations.
DoF can be queried for a given atom with `ctx.get_dof(atom_index)` and `ctx.get_dof_directional(atom_index)`,
or in bulk using `ctx.get_all_dof()` and `ctx.get_all_dof_directional()`.
The latter two functions can optionally take a parameter for a list of particular atoms to query,
`atoms=np.array(list_of_atoms, dtype=np.int64_t)`, and one for a buffer in which to store the output,
`buffer=np.zeros((len(list_of_atoms),))` or `buffer=np.zeros((len(list_of_atoms,3)))`.

Iterators over rigid and semi-rigid fragments can be obtained with `ctx.get_rigid_fragments()` and
`ctx.get_semirigid_fragments()`. For example:
```python
for rigid_body in ctx.get_rigid_fragments():
    for atom_index in rigid_body.atoms():
        do_something()
```
Note, fragment list iterators and fragment handles (elements returned by an iterator) become
invalid if the underlying data is modified, which can happen if `ctx.add_rigid_bond(...)`,
`ctx.build_rigid_fragment(...)`, or `ctx.finalise_fragments(...)` are called.


## MDAnalysis

To use dofulator with MDAnalysis, simply import the `MDADofulator` class,
construct an instance, and run the analysis.

For example:
```python
import MDAnalysis as mda
from dofulator import MDADofulator

# NOTE: Universe should include correct mass and position data,
#       Bond topology data is also useful to conveniently define constraints.
#       If the system is in a periodic volume, then periodic boundary information should also be present.
#       Velocity data is not required for DoF-only calculations, but is necessary for temperature.
u = mda.Universe('topol.tpr', 'traj.trr')

# Build list of rigid bodies (choosing residues with name 'H2O')
rigid_bodies = [r for r in u.residues if r.resname == 'H2O']

# List of rigid bonds. Should not contain any atoms in any of the rigid bodies.
# Elements of the list can be either `Bond` type or a 2-element iterable (e.g. a tuple) of `Atom`
# Here, choosing bond types 2, 12 and 20
rigid_bonds = [b for b in u.bonds if b.type in ['2', '12', 20']]

# Alternatively, to make all bonds to H rigid
rigid_bonds = [b for b in u.bonds if b[0].element == 'H' or b[1].element == 'H']

# Rigid angle list can be constructed similarly. Note that the A--B--C angle constraint
# is treated as three bond constraints between A--B, A--C and B--C, and again these should
# not intersect with atoms in rigid bodies.
# Here, H--X--H angles are being treated as constrained.
rigid_angles = [a for a in u.angles if a[0].element == 'H' and a[2].element == 'H']

# Lastly, rigid dihedrals work the same way as angles - by using bond constraints between
# A--B, A--C, A--D, B--C, B--D, and C--D in A--B--C--D
# If no rigid bodies/bonds/angles/dihedrals are present, `None` can be specified,
# which is the default.
rigid_dihedrals = None

d = MDADofulator(
    u.atoms,                            # Atoms to calculate DoF of
    rigid_bodies=rigid_bodies,          # Can alternatively pass in 'all' to treat all residues as rigid bodies
    rigid_bonds=rigid_bonds,            # As above, 'all' will treat all bonds as rigid
    rigid_angles=rigid_angles,          # As above, 'all' will treat all angles as rigid
    rigid_dihedrals=rigid_dihedrals,    # As above, 'all' will treat all dihedrals as rigid
    mode='atomic',          # (default) Only calculate total DoF of each atom
                            #   Alternatively, set mode='directional' to get x,y,z DoF
    verbose=True,           # (default), set `False` to disable progress bar (as for other MDA classes)
    null_space_thresh=None, # (default) - Use the C default value for null space threshold.
                            #   See above "Considerations for kinematic loops" for details of
                            #   what this does and typical values.
    use_pbc=True,           # (default) - Set `False` for non-periodic systems or to suppress warnings
                            #             about unknown periodic dimensions. If `True`, periodicity
                            #             is inferred from u.dimensions
    )

# Run the analysis
d.run()

# For mode='atomic':
# d.results now contains total DoF of each atom on each frame as a 2D numpy array.
# Index using `d.results[frame_id, atom_id]`, where `atom_id` is the index of the
# atom in the list of atoms passed in when `d` was constructed.

# For mode='directional':
# d.results is a 3D numpy array, indexed as `d.results[frame_id, atom_id, dir]`,
# where `dir` is 0, 1 or 2 for x, y or z DoF component.
```

To calculate local temperatures, a `LocalTemperature` analysis class is provided,
where each local temperature is defined by an
[MDAnalysis atom selection](https://docs.mdanalysis.org/stable/documentation_pages/selections.html).
For example, continuing from above:
```python
from dofulator import LocalTemperature
import numpy as np

# Build list of atom selections, each of which will get their temperature calculated on each frame.
# For example, using 10 equal size slabs in the z direction:
zrange = np.linspace(0, u.dimensions[2], 11) # 10 bins in the z direction
selections = [u.select_atoms(f'prop z >= {zrange[i]} and prop z < {zrange[i+1]}', updating=True)
        for i in range(len(zrange) - 1)],
# For spatial selections (or other selections which may change as the trajectory progresses),
# make sure to set `updating=True`!

t = LocalTemperature(
    selections,                 # List of atom selections to calculate temperature of.
    d,                          # Source of DoF - see below for details.
    store_dof_results=False,    # (default) Set True to also store the DoF values in d.results (slower)
    mode='atomic',              # (default) Calculate total temperature of each selection.
                                #   Set mode='directional' to calculate directional temperatures.
    verbose=True,               # (default) as above
    )

t.run()
```
Note, the atom list passed in to the `MDADofulator` class when it was constructed *MUST*
contain all atoms which may possibly be selected by any of the selections.

Care should be taken when using selections which may sometimes be empty (e.g. very small spatial bins).
Division by the zero DoF in those selections will result in `nan` or `inf` temperatures,
the handling of which should be carefully considered.

For some constraint topologies, the total DoF of a given atom does not vary greatly
as the molecule moves, and hence it can be useful for performance reasons to pre-calculate
DoF (e.g. from the equilibrium molecular geometry) and use a fixed value for `LocalTemperature`.
This approach works best when selections contain a relatively large number of the same molecule,
and are not biased towards a particular section of the molecule, so that any random fluctuations
will tend to cancel on average. However, care should be taken with inhomogeneous systems where
molecules may have different average geometries in different regions.

To use a fixed DoF value instead of an `MDADofulator` object, simply pass in a Numpy array containing
the DoF of each atom that might be included in the selections, indexed by `atom.ix`.
`d.result[0]` from the previous section would be suitable for this if the first frame contained molecules in their
equilibrium geometry.
