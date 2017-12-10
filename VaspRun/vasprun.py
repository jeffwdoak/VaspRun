#!/usr/bin/env python

"""
vasprun.py Jeff Doak jeff.w.doak@gmail.com

Tool to read in a vasprun.xml file and extract calculation information from it.

"""

import xml.etree.ElementTree as ET
import numpy as np

def search(node):
    """
    Function to list the branches of a node on the xml tree structure.
    Useful for searching stuff.

    """
    for i in node:
        print i.tag, i.attrib, i.text

class VaspRun(object):
    """
    Class to read in a vasprun.xml file and extract desired content.

    Attributes
    ----------
    tree : ElementTree
        The parsed vasprun.xml xml tree.
    root : Element
        The root element of the vasprun.xml tree.
    atom_names : list
        List of the names of every atom in the unitcell
    atom_types : list
        List of the types of every atom in the unitcell
    atom_type_names : list
        List of the name of each type of atom in the unitcell

    """

    def __init__(self, input_=None):
        if input_ is None:
            self.tree = ET.parse('vasprun.xml')
        else:
            self.tree = ET.parse(input_)
        self.root = self.tree.getroot()
        self.atom_names = [
            i[0].text for i in
            self.root.find('atominfo').findall('array')[0].find('set')
        ]
        self.atom_indices = [
            int(i[1].text)-1 for i in
            self.root.find('atominfo').findall('array')[0].find('set')
        ]
        self.atom_types = [
            int(i[0].text) for i in
            self.root.find('atominfo').findall('array')[1].find('set')
        ]
        self.atom_type_names = [
            i[1].text for i in
            self.root.find('atominfo').findall('array')[1].find('set')
        ]

    def recip_lat(self, flag="init"):
        """
        Read in reciprocal lattice vectors of initial calculation. Lattice
        vectors are missing a factor of 2*pi.

        """
        recip = np.zeros((3, 3))
        if flag == 'init':
            entry = self.root.findall('structure')[0].find('crystal').findall('varray')[1]
        else:
            entry = self.root.findall('structure')[1].find('crystal').findall('varray')[1]
        for i in range(3):
            temp = entry[i].text.split()
            for j in range(3):
                recip[i, j] = float(temp[j])
        return recip

    def read_dos(self, flag="total"):
        """
        Read in electronic DOS, either total or partial.

        """
        dos_input = self.root.find('calculation').find('dos')
        fermi_level = float(dos_input.find('i').text)
        # specific to total dos
        dos_array = dos_input.find(flag)[0][-1][-1]
        energy = np.zeros(len(dos_array))
        dos = np.zeros((len(dos_array), 2))
        for i, line in enumerate(dos_array):
            energy[i] = float(line.text.split()[0])
            dos[i, 0] = float(line.text.split()[1])
            dos[i, 1] = float(line.text.split()[2])
        return fermi_level, energy, dos

    def read_kpt_gen(self):
        """
        Read in data used to generate k-point mesh/path.

        """
        gen = self.root.find('kpoints').find('generation').findall('v')
        nkpts = int(self.root.find('kpoints').find('generation').find('i').text)
        kpath = np.zeros([len(gen), 3])
        for i in range(len(gen)):
            for j in range(3):
                kpath[i, j] = float(gen[i].text.split()[j])
        return kpath, nkpts

    def read_kpoints(self):
        """
        Read in k-points used in calculation.

        """
        kpoint_input = self.root.find('kpoints').findall('varray')[0]
        nkpts = len(kpoint_input)
        kpts = np.zeros((nkpts, 3))
        for i in range(nkpts):
            kpts[i, 0] = float(kpoint_input[i].text.split()[0])
            kpts[i, 1] = float(kpoint_input[i].text.split()[1])
            kpts[i, 2] = float(kpoint_input[i].text.split()[2])
        return kpts

    def read_kpt_centering(self):
        """
        Read in the generating scheme of the k-point mesh, either Gamma or
        Monkhorst-Pack.

        """
        return self.root.find('kpoints').find('generation').attrib['param']

    def read_kpt_weights(self):
        """
        Read in multiplicities of k-points used in calculation.

        """
        weight_input = self.root.find('kpoints').findall('varray')[1]
        nweights = len(weight_input)
        weights = np.zeros(nweights)
        for i in range(nweights):
            weights[i] = float(weight_input[i].text)
        return weights

    def read_kpt_projection(self):
        """
        Read in supercell k-points projected onto primitive-cell recpirocal
        lattice.

        """
        kpts = self.root.find('calculation').find('kprojected')
        # Read in k-point projection weights
        weights = kpts.find('array').find('set')
        nspins = len(weights)
        nsupkpts = len(weights[0])
        nbands = len(weights[0][0])
        nprimkpts = len(weights[0][0][0][0].text.split())
        projection = np.zeros((nspins, nsupkpts, nbands, nprimkpts, 3))
        # loop over spin, supercell k-points, and band index to read in the
        # projection weights
        for i in range(nspins):
            for j in range(nsupkpts):
                for k in range(nbands):
                    for m in range(nprimkpts):
                        projection[i, j, k, m, 2] = weights[i][j][k][0].text.split()[m]
        # Read in electron eigenvalues and occupations
        eigenvals = kpts.find('eigenvalues').find('array').find('set')
        for i in range(nspins):
            for j in range(nsupkpts):
                for k in range(nbands):
                    for m in range(nprimkpts):
                        projection[i, j, k, m, 0] = eigenvals[i][j][k].text.split()[0]
                        projection[i, j, k, m, 1] = eigenvals[i][j][k].text.split()[1]
        return projection

    def read_prim_kpts(self):
        """
        Read in kpoints of IRZ of primitive cell.

        """
        kpt_list = self.root.find('calculation').find('kpoints').findall('varray')[0]
        nkpts = len(kpt_list)
        kpts = np.zeros((nkpts, 3))
        for i in range(nkpts):
            for j in range(3):
                kpts[i, j] = float(kpt_list[i].text.split()[j])
        return kpts

    def read_eigenvals(self):
        """
        Read in electron eigenvalues.
        Dimensions of returned array are:
        # spins, # kpts, # bands, 2 (eigenvalue, occupation)

        """
        try:
            eigendata = self.root.find('calculation').find('projected').find('eigenvalues').find('array')[-1]
        except AttributeError:
            eigendata = self.root.find('calculation').find('eigenvalues').find('array')[-1]
        nspins = len(eigendata)
        nkpts = len(eigendata[0])
        nbands = len(eigendata[0][0])
        eigenvals = np.zeros((nspins, nkpts, nbands, 2))
        for i in range(nspins):  # loop over spins
            for j in range(nkpts):  # loop over k-points
                for k in range(nbands):  # loop over bands
                    eigenvals[i, j, k, 0] = float(eigendata[i][j][k].text.split()[0])
                    eigenvals[i, j, k, 1] = float(eigendata[i][j][k].text.split()[1])
        return eigenvals

    def read_eigenvec_projection(self):
        """
        Read in the projection of electron eigenvectors onto atomic orbitals.

        """
        projdata = self.root.find('calculation').find('projected').find('array')[-1]
        nspins = len(projdata)
        nkpts = len(projdata[0])
        nbands = len(projdata[0][0])
        nions = len(projdata[0][0][0])
        norbs = len(projdata[0][0][0][0].text.split())
        projection = np.zeros((nspins, nkpts, nbands, nions, norbs))
        for i in range(nspins):
            for j in range(nkpts):
                for k in range(nbands):
                    for l in range(nions):
                        for m in range(norbs):
                            projection[i, j, k, l, m] = float(projdata[i][j][k][l].text.split()[m])
        return projection
