import os
import sys
import pytest
import dofulator as dof
import MDAnalysis as mda
import numpy as np

DATADIR = os.path.join(os.path.dirname(__file__), 'data')
if not os.path.exists(DATADIR):
    print(f'Cannot find data directory: {DATADIR}', file=sys.stderr)
    exit(1)

TEMPDIR = os.path.join(DATADIR, 'temperature')
if not os.path.exists(TEMPDIR):
    print(f'Cannot find data directory: {TEMPDIR}', file=sys.stderr)
    exit(1)

"""
Testing of MDA interface
"""

@pytest.fixture
def bicyclo_conformations():
    """
    Trajectory with 1 molecule of 1-methoxybicyclo[2.2.1]heptane in different conformations.
    If treated as semi-rigid fragment, hits path for multiple loop closures, including
    overconstrained sections (rigid angles or dihedrals) and becoming fully rigid (rigid dihedrals).
    """
    return mda.Universe(
            os.path.join(DATADIR, '1-methoxybicyclo[2.2.1]heptane', 'geometry.gro'),
            os.path.join(DATADIR, '1-methoxybicyclo[2.2.1]heptane', 'conformations.trr'),
            to_guess=('masses', 'bonds', 'angles', 'dihedrals'),
        )

@pytest.fixture
def temperature_systems():
    """
    Systems with velocity data for testing local temperature calculation
    """
    return load_temperature_systems()

def load_temperature_systems():
    """
    Helper function to allow direct calling for data generation
    """
    return [
            (name, mda.Universe(
                os.path.join(DATADIR, name, 'topol.gro'),
                os.path.join(DATADIR, name, 'traj.trr'),
                to_guess=('masses', 'elements', 'bonds', 'angles', 'dihedrals'),
            ))
            for name in ('ethane', 'benzene')
        ]

class TestMDADofulator:
    def test_equivalent_inputs(self, bicyclo_conformations):
        """
        Test that rigid_angles is equivalent to also setting rigid_bonds for bonds in angles and
        rigid_dihedrals equivalent to also setting rigid_angles and rigid_bonds for those contained
        in the dihedrals.
        """
        bonds = [b for b in bicyclo_conformations.bonds if np.any([(b[0] in a and b[1] in a) for a in bicyclo_conformations.angles])]
        d_a = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles='all', use_pbc=False)
        d_ba = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds=bonds, rigid_angles='all', use_pbc=False)
        d_ba.run()
        d_a.run()
        np.testing.assert_array_almost_equal(d_a.results, d_ba.results, err_msg='Setting rigid_angles not equivalent to setting rigid_bonds + rigid_angles')

        bonds = [b for b in bicyclo_conformations.bonds if np.any([(b[0] in d and b[1] in d) for d in bicyclo_conformations.dihedrals])]
        angles = [a for a in bicyclo_conformations.angles if np.any([(a[0] in d and a[1] in d and a[2] in d) for d in bicyclo_conformations.dihedrals])]
        d_d = dof.MDADofulator(bicyclo_conformations.atoms, rigid_dihedrals='all', use_pbc=False)
        d_ad = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles=angles, rigid_dihedrals='all', use_pbc=False)
        d_bd = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds=bonds, rigid_dihedrals='all', use_pbc=False)
        d_bad = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds=bonds, rigid_angles=angles, rigid_dihedrals='all', use_pbc=False)

        d_d.run()
        d_ad.run()
        d_bd.run()
        d_bad.run()
        np.testing.assert_array_almost_equal(d_d.results, d_ad.results, err_msg='Setting rigid_dihedrals not equivalent to setting rigid_angles + rigid_dihedrals')
        np.testing.assert_array_almost_equal(d_d.results, d_bd.results, err_msg='Setting rigid_dihedrals not equivalent to setting rigid_bonds + rigid_dihedrals')
        np.testing.assert_array_almost_equal(d_d.results, d_bad.results, err_msg='Setting rigid_dihedrals not equivalent to setting rigid_bonds + rigid_angles + rigid_dihedrals')

    def test_two_loops(self, bicyclo_conformations):
        """
        1-methoxybicyclo[2.2.1]heptane with rigid bonds. Check total DoF is correct in all conformations.
        """
        d = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds='all', use_pbc=False)
        d.run()
        np.testing.assert_array_almost_equal(
                np.sum(d.results, axis=1),
                np.ones((d.results.shape[0],)) * \
                        (3*bicyclo_conformations.atoms.n_atoms - len(bicyclo_conformations.bonds)),
                err_msg='Incorrect total DoF with two kinematic loops'
            )

    def test_many_loops(self, bicyclo_conformations):
        """
        1-methoxybicyclo[2.2.1]heptane with rigid bonds and angles (many kinematic loops, overconstrained).
        Total DoF in all conformations should be 8:
         * 3 translational
         * 3 rotational
         * methoxy rotation around C--O (C in loop structure)
         * methyl rotation around O--C (C in methyl)
        """
        d = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles='all', use_pbc=False)
        d.run()
        np.testing.assert_array_almost_equal(
                np.sum(d.results, axis=1),
                np.ones((d.results.shape[0],)) * 8,
                err_msg='Incorrect total DoF with many kinematic loops + overconstrained sections'
            )

    def test_rigid_via_loops(self, bicyclo_conformations):
        """
        1-methoxybicyclo[2.2.1]heptane with rigid bonds, angles and dihedrals
        (many kinematic loops, overconstrained, equivalent to rigid body).
        """
        bicyclo_conformations.trajectory[1]
        d_semirigid = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles='all', rigid_dihedrals='all', use_pbc=False)
        d_rigid = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bodies='all', use_pbc=False)

        # Bond lengths/angles/dihedrals change each conformation, so check each frame individually.
        # Otherwise rigid case assumption of fixed geometry is incorrect.
        for _  in bicyclo_conformations.trajectory:
            d_semirigid.run_single_frame()
            d_rigid.run_single_frame()

            np.testing.assert_array_almost_equal(
                    d_semirigid.results, d_rigid.results,
                    err_msg="Atomic DoF from fully constrained semirigid fragment don't match those of rigid body",
                )

        d_semirigid.mode = 'directional'
        d_rigid.mode = 'directional'
        for _ in bicyclo_conformations.trajectory:
            d_semirigid.run_single_frame()
            d_rigid.run_single_frame()
            np.testing.assert_array_almost_equal(
                    d_semirigid.results, d_rigid.results,
                    err_msg="Directional DoF from fully constrained semirigid fragment don't match those of rigid body",
                )


class TestLocalTemperatureRegression:
    def test_temperature_configs(self, temperature_systems):
        """
        Check that calculated local temperatures match previously generated data files.
        """
        for (name, u) in temperature_systems:
            def check_results(t: dof.LocalTemperature, bonds, angles, bodies, mode):
                fname = get_temperature_fname(name, bonds, angles, bodies, mode)
                expected = np.fromfile(fname).reshape(t.results.shape)
                np.testing.assert_array_equal(t.results, expected)
            loop_temperature_constraint_combos(u, check_results)

def get_temperature_fname(name: None|str, bonds: None|str, angles: None|str, bodies: None|str, mode: None|str):
    """
    Get file of temperature data for given configuration.
    """
    fname = f'{name}.{mode}'
    if bonds is not None:
        fname += '.rigidbonds'
    if angles is not None:
        fname += '.rigidangles'
    if bodies is not None:
        fname += '.rigidbodies'
    fname += '.npy'
    return os.path.join(TEMPDIR, fname)


def loop_temperature_constraint_combos(u: mda.Universe, fn):
    """
    Calculate local temperatures for combinations of rigid bonds, angles and
    bodies in atomic and directional modes, and apply `fn` to each case.
    """
    d = dof.MDADofulator(u.atoms)
    selections = [u.select_atoms(f'element {e}') for e in np.unique(u.atoms.elements)]
    selections.insert(0, u.select_atoms('all'))
    for bonds in [None, 'all']:
        for angles in [None, 'all']:
            for bodies in [None, 'all']:
                if bodies is not None and (bonds is not None or angles is not None):
                    continue
                for mode in ['atomic', 'directional']:
                    d.set_rigid_bonds(bonds)
                    d.set_rigid_angles(angles)
                    d.set_rigid_bodies(bodies)
                    t = dof.LocalTemperature(selections, d, mode=mode)
                    t.boltz *= 100 # Account for incorrect velocity units from .trr
                    t.run()
                    fn(t, bonds, angles, bodies, mode)


if __name__ == '__main__':
    # Generate temperature data files for regression testing
    systems = load_temperature_systems()
    for (name, u) in systems:
        def write_results(t: dof.LocalTemperature, bonds, angles, bodies, mode):
            fname = get_temperature_fname(name, bonds, angles, bodies, mode)
            t.results.tofile(fname)
        loop_temperature_constraint_combos(u, write_results)

