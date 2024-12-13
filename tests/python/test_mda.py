import pytest
import dofulator as dof
import MDAnalysis as mda
import numpy as np

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
            'data/1-methoxybicyclo[2.2.1]heptane.gro',
            'data/1-methoxybicyclo[2.2.1]heptane.trr',
            to_guess=('masses', 'bonds', 'angles', 'dihedrals'),
        )

class TestMDADofulator:
    def test_equivalent_inputs(self, bicyclo_conformations):
        """
        Test that rigid_angles is equivalent to also setting rigid_bonds for bonds in angles and
        rigid_dihedrals equivalent to also setting rigid_angles and rigid_bonds for those contained
        in the dihedrals.
        """
        bonds = [b for b in bicyclo_conformations.bonds if np.any([(b[0] in a and b[1] in a) for a in bicyclo_conformations.angles])]
        d_a = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles='all')
        d_ba = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds=bonds, rigid_angles='all')
        d_ba.run()
        d_a.run()
        np.testing.assert_array_almost_equal(d_a.results, d_ba.results, err_msg='Setting rigid_angles not equivalent to setting rigid_bonds + rigid_angles')

        bonds = [b for b in bicyclo_conformations.bonds if np.any([(b[0] in d and b[1] in d) for d in bicyclo_conformations.dihedrals])]
        angles = [a for a in bicyclo_conformations.angles if np.any([(a[0] in d and a[1] in d and a[2] in d) for d in bicyclo_conformations.dihedrals])]
        d_d = dof.MDADofulator(bicyclo_conformations.atoms, rigid_dihedrals='all')
        d_ad = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles=angles, rigid_dihedrals='all')
        d_bd = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds=bonds, rigid_dihedrals='all')
        d_bad = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds=bonds, rigid_angles=angles, rigid_dihedrals='all')

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
        d = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bonds='all')
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
        d = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles='all')
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
        d_semirigid = dof.MDADofulator(bicyclo_conformations.atoms, rigid_angles='all', rigid_dihedrals='all')
        d_rigid = dof.MDADofulator(bicyclo_conformations.atoms, rigid_bodies='all')

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


class TestLocalTemperature:
    def test_trivial(self):
        return True;
