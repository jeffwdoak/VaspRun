#!/usr/bin/env python

# band_structure.py v1.0 5-08-2014 Jeff Doak jeff.w.doak@gmail.com

import sys
from vasprun import *
import numpy as np

vasprun = VaspRun()

recip = vasprun.recip_lat()
kpoints = vasprun.read_kpoints()
eigenvals = vasprun.read_eigenvals()
kpath,nkpts = vasprun.read_kpt_gen()

# Convert k-points to cartesian coordinates
kpoints = np.dot(kpoints,recip)
kpath = np.dot(kpath,recip)

# Generate k-point path, differentiating between breaks in the path
#path_shape = []
#for i in range(len(kpath)-1):
#    path_shape.append([])
#    if i > 0:
#        dk = path_shape[i-1][-1]
#    else:
#        dk = 0
#    kmag = np.linalg.norm(kpath[i+1]-kpath[i])
#    for j in range(nkpts):
#        path_shape[i].append(dk+kmag*float(j)/float(nkpts-1))
#path_shape = np.array(path_shape)

#tol = 1e-8
#path_shape2 = [[]]
#for i in range(len(kpath)-1):
#    if i > 0:
#        dk = path_shape2[-1][-1][-1]
#        if ((kpath[i]-kpath[i-1])<tol).all():
#        #if kpath[i] != kpath[i-1]:
#            path_shape2.append([])
#    else:
#        dk = 0
#    path_shape2[-1].append([])
#    kmag = np.linalg.norm(kpath[i+1]-kpath[i])
#    for j in range(nkpts):
#        path_shape2[-1][-1].append(dk+kmag*float(j)/float(nkpts-1))




# Construct list of band-structure k-point magnitudes
dk = [0]
k0 = kpoints[0]
i=1
for kpt in kpoints[1:]:
    if i == nkpts:
        i = 0
        dk.append(dk[-1])
    else:
        dk.append(dk[-1]+np.linalg.norm(kpt-k0))
    k0 = kpt
    i += 1

text = ""
for i in range(len(eigenvals)):  # loop over spins
    for j in range(len(eigenvals[0])):  # loop over kpoints
        text += str(dk[j])+" "
        for k in range(len(eigenvals[0,0])):  # loop over eigenvals
            text += str(eigenvals[i,j,k,0])+" "
        text += "\n"

print text
