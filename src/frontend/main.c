#include "dofulator.h"
#include "parse_molecule.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct IO {
  char* in;
  FILE* out;
} IO;

int parse_args(int argc, char* argv[], IO* io);

int main(int argc, char* argv[]) {
  IO io;
  int ret = parse_args(argc, argv, &io);
  if (io.in == NULL) goto cleanup;

  Molecule mol = parse_molecule_file(io.in);
  if (mol.atoms.n > 0) {
    Fragment frag;
    if (mol.bonds.n == 0) {
      frag = fragment_create_rigid(mol.atoms);
    } else {
      frag = fragment_create_semirigid(mol.atoms, mol.bonds);
    }
    if (frag) {
      fprintf(io.out, "# id\tlabel\tdof_total\tdof_x\t\tdof_y\t\tdof_z\n");
      for (unsigned i = 0; i < mol.atoms.n; ++i) {
        double dof = fragment_dof_atom(frag, i);
        double dof_x = fragment_dof_atom_dir(frag, i, (double[3]){1.0, 0.0, 0.0});
        double dof_y = fragment_dof_atom_dir(frag, i, (double[3]){0.0, 1.0, 0.0});
        double dof_z = fragment_dof_atom_dir(frag, i, (double[3]){0.0, 0.0, 1.0});

        fprintf(io.out, "%d\t%-7s\t%.6f\t%.6f\t%.6f\t%.6f\n", i+1, mol.atom_labels[i], dof, dof_x, dof_y, dof_z);
      }
    } else {
      fprintf(stderr, "Error processing fragment\n");
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

  *io = (IO){.in = NULL, .out = stdout};
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
