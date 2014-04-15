#!/usr/bin/env python

# projected_bands_from_xml.py v1.0 3-20-2014 Jeff Doak jeff.w.doak@gmail.com

import sys
import numpy as np
from vasprun import *

def atom_projection(proj,atom_indices):
    """
    Function to sum the projection of each {kpoint,band} onto each atomic
    orbital across all orbitals of each atom of the same type. Returns an array
    of dimensions [nkpts,nbands,ntypes] containing the projection of {kpt,band}
    i,j onto atom k (i in 1 ... nkpts, j in 1 ... nbands, k in 1 ... ntypes).
    """
    ntypes = len(set(atom_indices))
    ai = atom_indices
    nspins = len(proj)
    nkpts = len(proj[0])
    nbands = len(proj[0][0])
    nions = len(proj[0][0][0])
    #atom_proj = np.zeros((nkpts,nbands,nions))
    atom_proj = np.zeros((nkpts,nbands,ntypes))
    for s in range(nspins):
        for i in range(nkpts):
            for j in range(nbands):
                for k in range(nions):
                    atom_proj[i,j,ai[k]] += np.sum(proj[s,i,j,k])
    return atom_proj

def atom_orbital_projection(proj,atom_indices):
    """
    Function to sum the projection of each {kpoint,band} onto each atomic
    orbital across all atoms of each type, keeping each of the orbital
    contributions distinct. Returns an array of dimensions
    [nkpts,nbands,ntypes,norbs] containing the projection of {kpoint,band} i,j
    onto orbital l (el) of atom k (i in 1 ... nkpts, j in 1 ... nbands, k in
    1 ... ntypes, l [el] in 1 ... norbs).
    """
    ntypes = len(set(atom_indices))
    ai = atom_indices
    nspins = len(proj)
    nkpts = len(proj[0])
    nbands = len(proj[0][0])
    nions = len(proj[0][0][0])
    norbs = len(proj[0][0][0][0])
    if norbs == 1:
        # s
        l_list = [1]
    elif norbs == 4:
        # p
        l_list = [1,3]
    elif norbs == 9:
        # d
        l_list = [1,3,5]
    elif norbs == 16:
        # f
        l_list = [1,3,5,7]
    else:
        print "Unknown # of orbitals"
        sys.exit(1)
    atom_orb_proj = np.zeros((nkpts,nbands,ntypes,len(l_list)))
    for s in range(nspins):
        for i in range(nkpts):
            for j in range(nbands):
                for k in range(nions):
                    n = 0
                    for l in range(len(l_list)):
                        for m in range(l_list[l]):
                            atom_orb_proj[i,j,ai[k],l] += proj[s,i,j,k,n]
                            n += 1
    return atom_orb_proj

if __name__ == "__main__":
    # Read in vasprun.xml data
    vasprun = VaspRun()
    kpts = vasprun.read_kpoints()
    eigenvals = vasprun.read_eigenval()
    projs = vasprun.read_eigenvec_projection()

    if len(sys.argv) > 1 and str(sys.argv[1]) == "atom":
        atom_proj_data = atom_projection(projs,vasprun.atom_indices)
    else:
        atom_proj_data = atom_orbital_projection(projs,vasprun.atom_indices)

    ntypes = len(vasprun.atom_types)
    nkpts = len(kpts)
    nbands = len(eigenvals[0][0])
    nions = len(projs[0][0][0])
    norbs = len(projs[0][0][0][0])

    # Construct list of k-point magnitudes
    delta_kpts = np.zeros(nkpts)
    for i in range(1,nkpts):
        delta_kpts[i] = np.linalg.norm(kpts[i]-kpts[i-1]) + delta_kpts[i-1]

    # Create and print header
    if norbs == 1:
        l_list = ['s']
    elif norbs == 4:
        l_list = ['s','p']
    elif norbs == 9:
        l_list = ['s','p','d']
    elif norbs == 16:
        l_list = ['s','p','d','f']
    line = "d_kpt eigenval_(eV) "
    for atom in vasprun.atom_type_names:
        for i in l_list:
            line += atom+"_"+i+" "
    print line

    # Loop over k-points and bands, printing out the projection onto each atom type
    # and orbital
    for i in range(nkpts):
        for j in range(nbands):
            line = str(delta_kpts[i])+" "
            line += str(eigenvals[0,i,j,0])+" "
            line += " ".join([ str(x) for x in atom_proj_data[i,j].flatten() ])
            print line
    sys.exit()
