#include "dofulator.h"
#include "parse_molecule.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

typedef struct IO {
  char* in;
  FILE* out;
  bool all_rigid;
} IO;

int parse_args(int argc, char* argv[], IO* io);

int main(int argc, char* argv[]) {
  IO io;
  int ret = parse_args(argc, argv, &io);
  if (io.in == NULL) goto cleanup;

  Molecule mol = parse_molecule_file(io.in);
  if (mol.atoms.n > 0) {
    Dofulator ctx = dofulator_create(mol.atoms.n);
    if (mol.bonds.n == 0) {
      for (size_t i = 1; i < mol.atoms.n; ++i) {
        dofulator_build_rigid_fragment(ctx, (Bond){0, i});
      }
    } else {
      for (Bond* b = mol.bonds.bonds; b < mol.bonds.bonds + mol.bonds.n; ++b) {
        if (io.all_rigid) {
          dofulator_build_rigid_fragment(ctx, *b);
        } else {
          dofulator_add_rigid_bond(ctx, *b);
        }
      }
    }
    dofulator_finalise_fragments(ctx);
    dofulator_precalculate_rigid(ctx, mol.atoms.mass, mol.atoms.x);
    dofulator_calculate(ctx, mol.atoms.mass, mol.atoms.x);
    if (ctx) {
      fprintf(io.out, "# id\tlabel\tdof_total\tdof_x\t\tdof_y\t\tdof_z\n");
      for (AtomTag i = 0; i < mol.atoms.n; ++i) {
        double dof = dofulator_get_dof_atom(ctx, i);
        double dof_dir[3];
        dofulator_get_dof_atom_directional(ctx, i, dof_dir);

        fprintf(io.out, "%ld\t%-7s\t%.6f\t%.6f\t%.6f\t%.6f\n", i+1, mol.atom_labels[i], dof, dof_dir[0], dof_dir[1], dof_dir[2]);
      }
    } else {
      fprintf(stderr, "Error building dofulator context\n");
      ret = EXIT_FAILURE;
    }
    molecule_destroy(&mol);
  } else {
    ret = EXIT_FAILURE;
  }

cleanup:
  if (io.out != stdout) fclose(io.out);

  return ret;
}


int parse_args(int argc, char* argv[], IO* io) {
  const char* USAGE =
    "Calculates total and directional degrees of freedom"
    "per atom, given a (modified) .xyz file format.\n"
    "\n"
    "USAGE:\n"
    "  dof   INPUT_FILE [-o/--output OUTPUT_FILE] [-h/--help]]\n"
    "\n"
    "FORMAT:\n"
    "  Input file for rigid molecules should be formatted as below.\n"
    "  Atom labels are optional, and must not exceed 7 characters.\n"
    "  Additional columns/text to the right of the required ones\n"
    "  may be present, and will be ignored.\n"
    "\n"
    "  [num_atoms]\n"
    "  [comment line (may be empty)]\n"
    "  [atom_id] [atom label] [mass] [x] [y] [z]\n"
    "  [...]\n"
    "\n\n"
    "  For semi-rigid fragments with rigid bond lengths but flexible\n"
    "  angles, the input file should additionally contain:\n"
    "\n"
    "  [num_bonds]\n"
    "  [atom1_id] [atom2_id]\n"
    "  [...]\n"
    "\n";

  *io = (IO){.in = NULL, .out = stdout, .all_rigid = false};
  for (int iarg = 1; iarg < argc; ++iarg) {
    if (strcmp("-o", argv[iarg]) == 0 || strcmp("--output", argv[iarg]) == 0) {
      if (++iarg >= argc) {
        fprintf(stderr, "ERROR: Expected output file name after -o/--output flag\n");
        return EXIT_FAILURE;
      }
      if (io->out != stdout) {
        fprintf(stderr, "WARNING: Multiple output files specified! Defaulting to first.\n");
      } else {
        io->out = fopen(argv[iarg], "w");
      }

    } else if (strcmp("--all-rigid-mols", argv[iarg]) == 0) {
      io->all_rigid = true;

    } else if (strcmp("-h", argv[iarg]) == 0 || strcmp("--help", argv[iarg]) == 0) {
      fprintf(stdout, "%s", USAGE);
      if (io->out != stdout) {
        fclose(io->out);
      }
      io->in = NULL;
      return EXIT_SUCCESS;

    } else if (io->in == NULL) {
      io->in = argv[iarg];

    } else {
      fprintf(stderr, "ERROR: Unexpected argument '%s'\n", argv[iarg]);
      return EXIT_FAILURE;
    }
  }

  return EXIT_SUCCESS;
}
