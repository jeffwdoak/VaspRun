#!/usr/bin/env python

# band_structure.py v1.0 5-08-2014 Jeff Doak jeff.w.doak@gmail.com

import sys
from vasprun import *
import numpy as np

vasprun = VaspRun()

recip = vasprun.recip_lat()
kpoints = vasprun.read_kpoints()
eigenvals = vasprun.read_eigenval()

kpoints = np.dot(kpoints,recip)  # Convert k-points to cartesian coordinates

# Construct list of band-structure k-point magnitudes
dk = [0]
k0 = kpoints[0]
for kpt in kpoints[1:]:
    dk.append(dk[-1]+np.linalg.norm(kpt-k0))
    k0 = kpt

text = ""
for i in range(len(eigenvals)):  # loop over spins
    for j in range(len(eigenvals[0])):  # loop over kpoints
        text += str(dk[j])+" "
        for k in range(len(eigenvals[0,0])):  # loop over eigenvals
            text += str(eigenvals[i,j,k,0])+" "
        text += "\n"

print text
