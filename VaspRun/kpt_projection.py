#!/usr/bin/env python

# kpt_projection.py v1.2 4-01-2014 Jeff Doak jeff.w.doak@gmail.com

# script to read in k-point projection scheme data from one or more vasprun.xml
# files and assemble an unfoleded 'effective' band structure for the
# calculation(s) in an associated primtive cell.

# Steps in the algorithm:
# 1. Identify files for all calculations involved in supercell k-points.
# 2. Read in kpt-projection scheme for each calculation's vasprun.xml file.
# 3. Concatenate primtive cell k-points and spectral weights.
# 4. ???
# 5. Profit

from vasprun import *
from kpt_transform import *
import numpy as np
import sys

# Implementation


def find_projected_kpts(k0,kpts,lat,sym):
    """
    Search through the list of primitive-cell-projected kpoints, kpts, to find 
    all the primitive-cell-projected kpoints which are symmetrically equivalent
    to the primtive-cell kpoint k0, under symmetry operations sym corresponding
    to the reciprocal lattice vectors, lat. What exactly should this function
    return? We'll find out when I use it later. For now it will just return a
    list of indices to back out the corresponding kpoints and kpoint
    projections. Also, I will need to think some about how to best implement
    the dirs part of the function.
    """
    list_ = []
    symkpts = remove_dup(in_ws_cell(apply_sym_ops(sym,[k0]),lat))
    for i in range(len(kpts)):
        for j in range(len(symkpts)):
            if (kpts[i] == symkpts[j]).all():
                list_.append(i)
    return list_


def spectral_function(kpt,prim_kpts,projection,lat,sym):
    """
    Function to calculate the spectral distrubution function for the projection
    of a supercell bandstructure onto the primitive cell bandstructure at a
    given primitive cell k-point, kpt. Returns an array containing pairs of
    electron eigenvalue and spectral weight at that eigenvalue (only non-zero
    spectral weights are returned).
    """
    symkpts = remove_dup(in_ws_cell(apply_sym_ops(sym,[kpt]),lat))
    # shape of projection:
    # nspins, nsupkpts, nbands, nprimkpts, [eng,occ,P_km]

    # Spectral function for kpt k0 is an array of length nbands
    nbands = len(projection[0][0,0])
    #spectral = np.zeros((nbands,2))
    spectral = []
    nkpts = 0

    tol = 1e-12
    # Loop over each supercell calculation
    for i in range(len(prim_kpts)):
        # Loop over each primitive-cell kpt in supercell calc.
        for j in range(len(prim_kpts[i])):
            # Loop over each kpt symetrically equiv. to k0
            for k in range(len(symkpts)):
                # See if kpt prim_kpts[i][j] is equiv to k0
                if (abs(prim_kpts[i][j] - symkpts[k]) < tol).all():
                    nkpts += 1
                #if (prim_kpts[i][j] == symkpts[k]).all():
                    # If true, add supercell kpoint projection weights to
                    # spectral function of k0
                    # Loop over all spins
                    for s in range(len(projection[i])):
                        # Loop over all supercell kpoints
                        for K in range(len(projection[i][s])):
                            # Loop over all bands
                            for m in range(len(projection[i][s,K])):
                                if projection[i][s,K,m,j,2] > 0:
                                    # Check to see if eigenvalue is already in
                                    # spectral function. If so, add weight to
                                    # that array. Otherwise add new entry.
                                    e = 0
                                    flag = 0
                                    while e < len(spectral):
                                        if abs(spectral[e][0] - projection[i][s,K,m,j,0]) < tol:
                                            spectral[e][1] += projection[i][s,K,m,j,2]
                                            flag = 1
                                            break
                                        else:
                                            e += 1
                                    if flag == 0:
                                        spectral.append([projection[i][s,K,m,j,0],projection[i][s,K,m,j,2]])
    spectral = np.array(spectral)
    #print np.shape(spectral)
    #spectral[:,1] = spectral[:,1]/float(len(symkpts))
    spectral[:,1] = spectral[:,1]/float(nkpts)
    return spectral

# New approach

# 1. Read in POSCAR and vasprun.xml files for generating primtive cell
#    calculation
sym_ops = read_sym_ops(str(sys.argv[1]))
genrun = VaspRun()
gen_kpts = genrun.read_kpoints()
gen_kpts = remove_dup(gen_kpts)
primcell = UnitCell('POSCAR.prim')
recip = primcell.recip_lat()

# 2. Read in vasprun.xml files for all supercell calculations
dirs = [ os.getcwd()+os.sep+i for i in os.listdir('.') if os.path.isdir(i) ]
sup_kpts = []
prim_kpts = []
projection = []
for i in dirs:
    vasprun = VaspRun(i+os.sep+'vasprun.xml')
    sup_kpts.append(vasprun.read_kpoints())
    prim_kpts.append(vasprun.read_prim_kpts())
    projection.append(vasprun.read_kpt_projection())

# Construct list of band-structure k-point magnitudes
dk = [0]
k0 = gen_kpts[0]
for kpt in gen_kpts[1:]:
    dk.append(dk[-1]+np.linalg.norm(kpt-k0))
    k0 = kpt

# Construct array of spectral weights for each generating primitive-cell kpoint
spec = []
for i,kpt in enumerate(gen_kpts):
    temp_spec = spectral_function(kpt,prim_kpts,projection,recip,sym_ops)
    if len(temp_spec) > 0:
        dk_array = np.ones(len(temp_spec))*dk[i]
        #temp_spec = np.insert(temp_spec,0,dk_array,1)
        spec.append(np.insert(temp_spec,0,dk_array,1))
spec = np.array(spec)

# Print out array of band-structure k-point magnitudes and spectral functions.
for line in spec:
    for i in line:
        print i[0],i[1],i[2]

sys.exit()
