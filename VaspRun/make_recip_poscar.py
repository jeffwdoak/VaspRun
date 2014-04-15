#!/usr/bin/env python

# make_recip_poscar.py v1.0 2-17-2014 Jeff Doak jeff.w.doak@gmail.com

import sys
import numpy as np
from unitcell import *
from vasprun import *

# Read in real-space unit cell, and calculate reciprocal lattice from it
# Using convention of vasp to not include factor of 2*pi in the reciprocal
# lattice vectors. Not sure how including or not including the factor of 2*pi
# will affect my results. I imagine it won't matter for the determiniation of
# symmetry elements.
real_space_cell = UnitCell(str(sys.argv[1]))
recip_lat = real_space_cell.recip_lat()

supercell = UnitCell(real_space_cell)
supercell.simple_supercell(3,3,3)


# Read in vasprun.xml file to get the list of k-points used in the primitive
# cell calculation.
vasprun = VaspRun()
kpts = vasprun.read_kpoints()

# Create new unitcell instance with the reciprocal lattice as the unit cell
# vectors, and kpoints as the atomic positions. This should create a
# 'poscar-like-file' which encapsulates the reciprocal lattice. I will then use
# this file as an input to mint, to apply symmetry operations.
names = []
for i in range(len(kpts)):
    names.append(None)
recip_space_cell = UnitCell()
recip_space_cell.name = "reciprocal_unit_cell_w_kpts"
recip_space_cell.cell_vec = recip_lat
recip_space_cell.num_atoms = len(kpts)
recip_space_cell.num_atom_types = 1
recip_space_cell.atom_types = [len(kpts)]
recip_space_cell._atom_type_names = ['kpt']
recip_space_cell._atom_names = names
recip_space_cell._convention = "Direct"
recip_space_cell.atom_positions = kpts

recip_cell_str = recip_space_cell.output_vasp()
print recip_cell_str
sys.exit()
