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
  Dofulator ctx = NULL;
  DofulatorResult err;

  int ret = parse_args(argc, argv, &io);
  if (io.in == NULL) goto cleanup;
  Molecule mol = parse_molecule_file(io.in);

  if (mol.atoms.n > 0) {
    ctx = dofulator_create(mol.atoms.n);
    if (mol.bonds.n == 0) {
      for (size_t i = 1; i < mol.atoms.n; ++i) {
        err = dofulator_build_rigid_fragment(ctx, (Bond){0, i});
        if (err) {
          fprintf(stderr, "Fragment building failed with error code %d\n", err);
          ret = err;
          goto cleanup;
        }
      }
    } else {
      for (Bond* b = mol.bonds.bonds; b < mol.bonds.bonds + mol.bonds.n; ++b) {
        if (io.all_rigid) {
          err = dofulator_build_rigid_fragment(ctx, *b);
        } else {
          err = dofulator_add_rigid_bond(ctx, *b);
        }
        if (err) {
          fprintf(stderr, "Fragment building failed with error code %d\n", err);
          ret = err;
          goto cleanup;
        }
      }
    }
    err = dofulator_finalise_fragments(ctx);
    if (err) {
      fprintf(stderr, "Error finalising fragments: code %d\n", err);
      ret = err;
      goto cleanup;
    }
    err = dofulator_precalculate_rigid(ctx, mol.atoms.mass, mol.atoms.x);
    if (err) {
      fprintf(stderr, "Error precalculating rigid body DoF: code %d\n", err);
      ret = err;
      goto cleanup;
    }
    err = dofulator_calculate(ctx, mol.atoms.mass, mol.atoms.x);
    if (err) {
      fprintf(stderr, "Error calculating DoF: code %d\n", err);
      ret = err;
      goto cleanup;
    }
    if (ctx) {
      fprintf(io.out, "# id\tlabel\tdof_total\tdof_x\t\tdof_y\t\tdof_z\n");
      for (AtomTag i = 0; i < mol.atoms.n; ++i) {
        double dof = dofulator_get_dof_atom(ctx, i);
        double dof_dir[3];
        dofulator_get_dof_atom_directional(ctx, i, dof_dir);

        fprintf(io.out, "%zu\t%-7s\t%.6f\t%.6f\t%.6f\t%.6f\n", i+1, mol.atom_labels[i], dof, dof_dir[0], dof_dir[1], dof_dir[2]);
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
  dofulator_destroy(&ctx);

  return ret;
}


int parse_args(int argc, char* argv[], IO* io) {
  const char* USAGE =
    "USAGE:\n"
    "  %s INPUT_FILE [-o/--output OUTPUT_FILE] [--all-rigid-mols] [-h/--help]]\n"
    "\n"
    "Calculates total and directional degrees of freedom"
    "per atom, given a (modified) .xyz file.\n"
    "\n"
    " ARG              DESCRIPTION\n"
    "  -o/--output      Name of output file. Will contain the total, x, y,\n"
    "                   and z DoF of each atom. Prints to stdout if not set.\n"
    "\n"
    "  --all-rigid-mols Treat all molecules (atoms connected by at least 1 bond)\n"
    "                   as rigid bodies.\n"
    "\n"
    "FILE FORMAT:\n"
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
    "\n"
    "  If the bonds section is not present, it is assumed that all\n"
    "  atoms form a single rigid body (otherwise each would have 3 DoF).\n"
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
      fprintf(stdout, USAGE, argv[0]);
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
