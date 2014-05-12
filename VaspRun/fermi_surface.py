#!/usr/bin/env python

# fermi_surface.py v1.0 5-01-2014 Jeff Doak jeff.w.doak@gmail.com

import sys
import numpy as np
from vasprun import *
from kpt_transform import *

# Read in data from vasprun.xml file
vasprun = VaspRun(str(sys.argv[1]))
eigenvals = vasprun.read_eigenvals()
kpoints = vasprun.read_kpoints()
fermi,energy,dos = vasprun.read_dos()
recip = vasprun.recip_lat()  # Missing factor of 2*pi. Not sure if that matters.

# Read in symmetry operations from file
real_ops = read_sym_ops(str(sys.argv[2]))
recip_ops = [ i.T for i in real_ops ]

# Check that origin of k-point mesh is Gamma
origin = vasprun.read_kpt_centering()
if origin != 'Gamma':
    print "K-point mesh not centered at Gamma! Cannot use this mesh with XCrySDen."
    print "Program will now exit!"
    sys.exit(1)

# Calculate # of kpoints along a* b* and c*
#kx = sorted(set(kpoints[:,0]))
#ky = sorted(set(kpoints[:,1]))
#kz = sorted(set(kpoints[:,2]))

#nkx = len(kx)
#nky = len(ky)
#nkz = len(kz)


# Calculate k-point mesh across full reciprocal unit cell
nkx,nky,nkz = [ int(i) for i in vasprun.root.find('kpoints').find('generation')[0].text.split() ]

def kpt(ikx,iky,ikz):
    """
    Function to calculate the fractional coordinates of a k-point corresponding
    to a set of k-point integers ranging from 0 to nkj-1, where j = x,y,z.
    """
    kx = ikx*1./float(nkx)
    ky = iky*1./float(nky)
    kz = ikz*1./float(nkz)
    return np.array((kx,ky,kz))

def irred_kpt(ikx,iky,ikz):
    """
    Function to determine the k-point in the irreducible BZ that corresponds to a
    given set of k-point integers ranging from 0 to nkj-1 where j = x,y,z.
    """
    k0 = kpt(ikx,iky,ikz)
    sym_equiv = remove_dup(apply_sym_ops(recip_ops,in_ws_cell([k0],recip)))
    sym_equiv -= sym_equiv.round()
    return sym_equiv
    #reducible = incell(sym_equiv)
    #norms = [ np.linalg.norm(i) for i in reducible ]
    #return reducible[np.argmin(norms)]

def kptnum(ikx,iky,ikz):
    """
    Function to map a set of kpt integers spanning the full reciprocal unit cell
    onto the index of the corresponding k-point in the irreducible BZ in the
    list kpoints, where eigenvalues are defined.
    """
    tol = 1e-08
    irred = irred_kpt(ikx,iky,ikz)
    #kptnum = [ i for i in range(len(kpoints)) for j in irred if (np.abs(kpoints[i]-j)<tol).all() ]
    for i in range(len(kpoints)):
        for j in irred:
            if (np.abs(kpoints[i]-j)<tol).all():
                break
    return i
    #return kptnum[0]
    #kptnum = [ i for i in range(len(kpoints)) if (np.abs(kpoints[i]-irred)<tol).all() ]
    #return kptnum


# Lets do the mapping the other direction
def meshpoint(k0):
    tol = 1e-8
    for ikx in range(nkx):
        if (np.abs(k0[0]-kpt(ikx,0,0)[0])<tol):
            break
    for iky in range(nky):
        if (np.abs(k0[1]-kpt(0,iky,0)[1])<tol):
            break
    for ikz in range(nkz):
        if (np.abs(k0[2]-kpt(0,0,ikz)[2])<tol):
            break
    return ikx,iky,ikz



kptmap = dict()
for i,k0 in enumerate(kpoints):
    sym_equiv = remove_dup(incell(apply_sym_ops(recip_ops,[k0])))
    for ksym in sym_equiv:
        ikx,iky,ikz = meshpoint(ksym)
        kptmap[(ikx,iky,ikz)] = i

# Create dictionary mapping kpt integers in the reciprocal primtive cell to the
# kpt # in the irreducible BZ.
#kptmap = dict()
#for i in range(nkx):
#    for j in range(nky):
#        for k in range(nkz):
#            kptmap[(i,j,k)] = kptnum(i,j,k)

# Reshape eigenvalues to match k-point mesh dimenions
eigenvals = eigenvals[0,:,:,0]
nbands = len(eigenvals[0])

# Write header for band structure file
text = "BEGIN_INFO\n"
text += "  Fermi Energy: "+str(fermi)+"\n"
text += "END_INFO\n"
text += "\n"
text += "BEGIN_BLOCK_BANDGRID_3D\n"
text += "  vasp_band_structure_data\n"
text += "  BEGIN_BANDGRID_3D\n"
text += "    "+str(nbands)+"\n"
#text += "    "+str(nkx)+" "+str(nky)+" "+str(nkz)+"\n"  # off-by one!
text += "    "+str(nkx+1)+" "+str(nky+1)+" "+str(nkz+1)+"\n"
text += "    0.0 0.0 0.0\n"  # origin of k-point mesh. Required to be Gamma.
# Write reciprocal lattice vectors to file
for i in range(3):
    text += "    "
    for j in range(3):
        text += str(recip[i,j])+" "
    text += "\n"
# Write bands to file
kxarray = range(nkx+1)
kxarray[-1] = 0
kyarray = range(nky+1)
kyarray[-1] = 0
kzarray = range(nkz+1)
kzarray[-1] = 0
for i in range(nbands):
    text += "  BAND: "+str(i+1)+"\n"
#    for ikx in range(nkx+1):
    for ikx in kxarray:
#        for iky in range(nky+1):
        for iky in kyarray:
            text += "      "
#            for ikz in range(nkz+1):
            for ikz in kzarray:
                text += str(eigenvals[kptmap[(ikx,iky,ikz)],i])+" "
            text+= "\n"
        text += "\n"
# old write bands to file 5/5/2014
#for i in range(nbands):
#    text += "  BAND: "+str(i+1)+"\n"
#    for kx in range(nkx):
#        for ky in range(nky):
#            text += "      "
#            for kz in range(nkz):
#                text += str(eigenvals[kx,ky,kz,i])+" "
#            text+= "\n"
#        text += "\n"
# Write footer for band structure file
text += "  END_BANDGRID_3D\n"
text += "END_BLOCK_BANDGRID_3D\n"

# Output text to file.
outfile = open('bandstructure.bxsf','w')
outfile.write(text)
outfile.close()

sys.exit()
