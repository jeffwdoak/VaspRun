#!/usr/bin/env python

# kpt_transform.py v1.0 3-22-2014 Jeff Doak jeff.w.doak@gmail.com

import sys
import numpy as np
from unitcell import *
from vasprun import *

def read_sym_ops(sym_name):
    """
    Reads in symmetry operations from a mint matrix-formatted output. Input
    variable sym_name should be a string containing the name of the file.
    """
    sym_file = open(sym_name,'r')
    nops = int(sym_file.readline().split()[-1])
    sym_ops = np.zeros((nops,3,3))
    i = 0
    j = 0
    while i < nops:
        line = sym_file.readline()
        if j != 3 and j != 4:
            line = line.split()
            for k in range(3):
                sym_ops[i,j,k] = float(line[k])
        i = i + (j+1)/5
        j = (j+1)%5
    sym_file.close()
    return sym_ops

def apply_sym_ops(ops,kpts):
    """
    Applies a set of symmetry operations to a set of k-points.
    """
    full_kpts = []
    for i in range(len(ops)):
        for j in range(len(kpts)):
            full_kpts.append(np.dot(ops[i],kpts[j]))
    full_kpts = np.array(full_kpts)
    return full_kpts

def trans_mat(prim,sup):
    """
    Function to calculate the transformation matrix that changes a basis from
    the primitive cell coordinates to the supercell coordinates. Arrays prim and
    sup are assumed to be 3x3 matricies containing the primitive and supercell
    lattice vectors as row vectors.
    """
    return np.dot(sup,np.linalg.inv(prim))

def incell(kpts,nmax=3):
    """
    Moves all kpts that lie outside of a reciprocal unit cell into the cell.
    Assumes that the array kpts is in fractional coordinates of some reciprocal
    lattice vectors.
    """
    tkpts = np.copy(kpts)
    i = nmax
    while i > 0:
        for j in range(len(tkpts)):
            for k in range(3):
                if tkpts[j,k] >= i:
                    tkpts[j,k] = tkpts[j,k] - i
                elif tkpts[j,k] < 1-i:
                    tkpts[j,k] = tkpts[j,k] + i
        i = i - 1
    return tkpts

def in_ws_cell(kpts,cell):
    """
    Moves all kpts that lie outside of the Wigner-Seitz cell of a reciprocal
    lattice defined by cell into the cell. Assumes that the array kpts are in
    fractional coordinates of the reciprocal lattice vectors in cell.
    """
    new_kpts = []
    met = cell.dot(cell.T)  # metrical matrix of lattice vectors
    for kpt in kpts:
        kpt -= np.round(kpt)
        dist = kpt.dot(cell)
        for j in range(3):
            dist -= np.round(dist.dot(cell[j])/met[j,j])*cell[j]
        new_kpts.append(dist.dot(np.linalg.inv(cell)))
    return np.array(new_kpts)

def kpt_mult(kpts,sym):
    """
    Determine the # of symmetrically equivalent k-points to each k-point in
    kpts, by applying the symmetry operations sym.
    """
    mult = np.zeros(len(kpts))
    for i in range(len(kpts)):
        tot_kpts = np.zeros([len(sym),3,3])
        for j in range(len(sym)):
            tot_kpts[j] = np.dot(sym[j],kpts[i])
        #for j in range(len(

def remove_dup(kpts):
    """
    Check list of kpts for duplicates, and return a list of kpts with all
    duplicates removed.
    """
    uniq_list = []
    tol = 1e-4
    for i in range(len(kpts)):
        flag = 0
        for j in range(i+1,len(kpts)):
            dk0 = np.abs(kpts[i,0] - kpts[j,0])
            dk1 = np.abs(kpts[i,1] - kpts[j,1])
            dk2 = np.abs(kpts[i,2] - kpts[j,2])
            #if dk0 > tol or dk1 > tol or dk2 > tol:
            if dk0 < tol and dk1 < tol and dk2 < tol:
                flag = 1
                break
        if flag == 0:
            uniq_list.append(kpts[i])
    uniq_list = np.array(uniq_list)
    return uniq_list

def printkpts(kpts):
    """
    Print a list of k-pts formatted as atomic positions in a cif file.
    """
    output = ""
    for i in range(len(kpts)):
        output += "K"+str(i+1)+" K "+str(kpts[i,0])+" "+str(kpts[i,1])+" "+str(kpts[i,2])+"\n"
    return output

def write_kpts_vasp(kpts):
    """
    Write a list of kpoints to a vasp-formatted KPOINTS file.
    """
    output = "Supercell_k-points_from_primitive_cell\n"
    output += str(len(kpts))+"\n"
    output += "reciprocal\n"
    for i in range(len(kpts)):
        for j in range(3):
            output += str(kpts[i,j])+" "
        output += "1.0\n"
    return output

if __name__ == '__main__':
    # Read in POSCARs of primitive- and super-cells
    prim  = UnitCell(str(sys.argv[1]))
    sup = UnitCell(str(sys.argv[2]))
    # Read in symmetry operations of primitive cell? Should it be supercell instead?
    sym_ops = read_sym_ops(str(sys.argv[3]))
    # set lattice constant of primitive and supercells to 1.0
    prim.set_scale()
    sup.set_scale()
    # calculate reciprocal lattice vectors of primitive and supercells
    prim_recip = prim.recip_lat()
    sup_recip = sup.recip_lat()
    # Find transformation matrix from primitve cell to supercell
    M = trans_mat(prim.cell_vec,sup.cell_vec)
    # Read in k-points of primitive-cell
    vasprun = VaspRun()
    prim_kpts = vasprun.read_kpoints()
    if str(sys.argv[-1]) != 'star':
        # 'Algorithm' to create set of k-points for a no-symmetry supercell calculation
        # that will allow for a statistical un-folding of the primitive BZ later on.

        # First, find all the k-points which are symmetrically equivalent to the
        # original set of k-points.
        # Note that the symmetry operations, sym_ops, need to be the symmetry operations in
        # the primitive reciprocal cell basis (e.g. bcc symmetry operations for fcc pc).
        all_prim_kpts = apply_sym_ops(sym_ops,prim_kpts)
        # Map primitive cell k-points into supercell BZ.
        all_sup_kpts = np.dot(M,all_prim_kpts.transpose()).transpose()
        # Fold supercell k-points into the supercell BZ.
        folded_sup_kpts = in_ws_cell(all_sup_kpts,sup_recip)
        # Remove duplicate (folded) k-points from the set of supercell k-points.
        uniq_folded_sup_kpts = remove_dup(folded_sup_kpts)
        # Write k-points to a vasp-formatted file
        nlines = int(sys.argv[-1])
        nfiles = len(uniq_folded_sup_kpts)/nlines
        rem = len(uniq_folded_sup_kpts)%nlines
        n = 0
        for i in range(nfiles):
            sup_kpt_file = open('KPOINTS_'+str(i),'w')
            output =  write_kpts_vasp(uniq_folded_sup_kpts[n:n+nlines])
            sup_kpt_file.write(output)
            sup_kpt_file.close()
            n += nlines
        if rem != 0:
            sup_kpt_file = open('KPOINTS_'+str(i+1),'w')
            output =  write_kpts_vasp(uniq_folded_sup_kpts[n:n+rem])
            sup_kpt_file.write(output)
            sup_kpt_file.close()

        sys.exit()
        # Check how things went.
        print "primitive k-pts"
        print len(all_prim_kpts)
        print
        print "supercell k-pts"
        print len(all_sup_kpts)
        print len(folded_sup_kpts)
        print
        #print folded_sup_kpts
        print len(uniq_folded_sup_kpts)
        print 
        print uniq_folded_sup_kpts
        sys.exit()
    else:
        # Calculate the star of each kpoint
        n = 0
        for kpt in prim_kpts:
            prim_star = apply_sym_ops(sym_ops,[kpt])
            sup_star = np.dot(M,prim_star.transpose()).transpose()
            folded_sup_star = in_ws_cell(sup_star,sup_recip)
            uniq_folded_sup_star = remove_dup(folded_sup_star)
            sup_kpt_file = open('KPOINTS_'+str(n),'w')
            output = write_kpts_vasp(uniq_folded_sup_star)
            sup_kpt_file.write(output)
            sup_kpt_file.close()
            n += 1
        sys.exit()




