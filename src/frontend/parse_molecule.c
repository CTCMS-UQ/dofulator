#include "parse_molecule.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

int skip_to_eol(FILE* fid) {
  for (int c; (c = getc(fid)) != '\n'; ) {
    if (c == EOF) return c;
  }
  return 0;
}

Molecule parse_molecule_file(char* fname) {
  Molecule out = (Molecule){.atoms = {.n = 0}, .bonds = {.n = 0}};

  // Open in binary mode to avoid fgetpos/ftell bug in MSVC standard library:
  // https://github.com/microsoft/STL/issues/1784
  FILE* fid = fopen(fname, "rb");
  if (!fid) {
    fprintf(stderr, "Error opening file:\n%s\n", fname);
    return out;
  }

  // Get number of atoms
  int n;
  if (fscanf(fid, "%d", &n) != 1) {
    fprintf(stderr, "Expected integer number of atoms as first line.\n");
    goto error;
  }
  if (n <= 0) {
    fprintf(stderr, "Expected a positive, non-zero number of atoms. Got: %d\n", n);
    goto error;
  }
  out.atoms.n = n;

  // Skip to end of atom number line and then end of description line
  if (skip_to_eol(fid) == EOF || skip_to_eol(fid) == EOF) {
    if (feof(fid)) {
      fprintf(stderr, "Unexpected EOF. No atoms present!\n");
      goto error;
    } else {
      fprintf(stderr, "Unexpected error skipping description line\n");
      goto error;
    }
  }

  // Allocate atom list and parse positions/masses
  out.atoms = atom_list_create(n);
  if (out.atoms.n == 0) {
    fprintf(stderr, "Error allocating atom list\n");
    goto error;
  }
  out.atom_labels = malloc(sizeof(*out.atom_labels) * n);
  int idx;
  char label[8] = "";
  double mass, x, y, z;
  unsigned idx_tot = 0;
  for (unsigned i = 0; i < out.atoms.n; ++i) {
    // Try without atom label first, or roll back to include one
    fpos_t fptr;
    fgetpos(fid, &fptr);
    if (fscanf(fid, "%d %lf %lf %lf %lf", &idx, &mass, &x, &y, &z) != 5) {
      if (fsetpos(fid, &fptr)) {
        fprintf(stderr, "Error retrying line with label");
        goto error;
      }
      if (fscanf(fid, "%d %7s %lf %lf %lf %lf", &idx, label, &mass, &x, &y, &z) != 6) {
        if (feof(fid)) {
          fprintf(stderr, "Unexpected end of file. Expected %d atoms, got %d.\n", n, i);
        } else {
          fprintf(stderr, "Invalid atom definition on line %d\n", i+3);
        }
        goto error;
      }
    } else {
      memset(label, 0, sizeof(label));
    }
    if (idx > n || idx == 0) {
      fprintf(fid, "Atom index out of range\n");
      goto error;
    }
    idx_tot += idx;
    --idx;  // Switch from 1-based to 0-based indexing
    out.atoms.x[idx][0] = x;
    out.atoms.x[idx][1] = y;
    out.atoms.x[idx][2] = z;
    out.atoms.mass[idx] = mass;
    memcpy(out.atom_labels[idx], label, sizeof(label));

    int eof_check = skip_to_eol(fid);
    if (i+1 < out.atoms.n && eof_check == EOF) {
      if (feof(fid)) {
        fprintf(stderr, "Unexpected end of file. Expected %d atoms, got %d.\n", n, i+1);
      } else {
        fprintf(stderr, "Unexpected error skipping to end of line %d.\n", i+3);
      }
      goto error;
    }
  }
  if (idx_tot != out.atoms.n*(out.atoms.n+1)/2) {
    fprintf(stderr, "Atom indices not unique\n");
    goto error;
  }

  // Look for bond list
  if (fscanf(fid, "%d", &n) != 1) {
    if (!feof(fid)) {
      fprintf(stderr, "Expected integer number of bonds or end of file.\n");
    }
    // If EOF, no bonds listed so just have rigid body
    goto done;
  }
  // Skip to and of bond number line
  if (skip_to_eol(fid) == EOF) {
    if (feof(fid)) {
      fprintf(stderr, "Unexpected EOF. No bonds present!\n");
      goto error;
    } else {
      fprintf(stderr, "Unexpected error reading bonds\n");
      goto error;
    }
  }

  // Allocate bond list and parse
  out.bonds = bond_list_create(n);
  if (out.bonds.n == 0) {
    fprintf(stderr, "Error allocating bond list\n");
    goto error;
  }
  int i, j;
  for (size_t b = 0; b < out.bonds.n; ++b) {
    if (fscanf(fid, "%d%d", &i, &j) != 2) {
      if (feof(fid)) {
        fprintf(stderr, "Unexpected end of file. Expected %zu bonds, got %zu.\n", out.bonds.n, b);
      } else {
        fprintf(stderr, "Expected bond as atom index pair. Invalid format for bond %zu\n", b+1);
      }
      goto error;
    }
    if (i > (int)out.atoms.n || i == 0 || j > (int)out.atoms.n || j == 0) {
      fprintf(stderr, "Bond atom index out of range\n");
      goto error;
    }
    // Switch from 1-based to 0-based indexing
    out.bonds.bonds[b] = (Bond){.ai = i-1, .aj = j-1};
    int eof_check = skip_to_eol(fid);
    if (b+1 < out.bonds.n && eof_check == EOF) {
      if (feof(fid)) {
        fprintf(stderr, "Unexpected end of file. Expected %zu bonds, got %zu.\n", out.bonds.n, b+1);
      } else {
        fprintf(stderr, "Unexpected error skipping to end of bond line.\n");
      }
      goto error;
    }
  }

  goto done;
error:
  molecule_destroy(&out);
done:
  fclose(fid);
  return out;
}


void molecule_destroy(Molecule* mol) {
  atom_list_destroy(&mol->atoms);
  bond_list_destroy(&mol->bonds);
  if (mol->atom_labels) {
    free(mol->atom_labels);
    mol->atom_labels = NULL;
  }
}
