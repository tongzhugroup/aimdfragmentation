#!/usr/bin/env python3
"""Fragmentation for AIMD"""

__author__ = "Jinzhe Zeng"
__email__ = "jzzeng@stu.ecnu.edu.cn"
__update__ = "2018-12-14"


import sys
import os
import time
from collections import Counter
from multiprocessing import Pool, cpu_count
import subprocess as sp
from ase.geometry import get_distances
from ase.io import read as readxyz
from ase import Atoms
import numpy as np
import openbabel
from GaussianRunner import GaussianRunner, GaussianAnalyst


class AIMDFragmentation(object):
    def __init__(self, nproc_sum=None, nproc=4, cutoff=3.5, xyzfilename="comb.xyz", pdbfilename="comb.pdb", qmmethod="mn15", qmbasis="6-31g(d)", addkw="", qmmem="400MW", atombondnumber=None, logfile=None, outputfile="force.dat", unit=1, pbc=False, cell=[0, 0, 0], gaussian_dir="gaussian_files", command="g16", gaussiancommand=None, jobfile="gaussianjobs", onebodykeyword="scf=(xqc,MaxConventionalCycles=256)", twobodykeyword="guess=mix scf=(maxcyc=256)", kbodyfile="kbforce.dat"):
        self.nproc_sum = nproc_sum if nproc_sum else cpu_count()
        self.nproc = nproc
        self.cutoff = cutoff
        self.xyzfilename = xyzfilename
        self.qmmethod = qmmethod
        self.qmbasis = qmbasis
        self.addkw = addkw
        self.qmmem = qmmem
        self.outputfile = outputfile
        self.unit = unit
        self.pbc = pbc
        self.cell = cell
        self._atomid = {}
        self.jobs = []
        self.gaussian_dir = gaussian_dir
        self.errorfiles = []
        self.command = command
        self.gaussiancommand = gaussiancommand
        self.jobfile = jobfile
        self.onebodykeyword = onebodykeyword
        self.twobodykeyword = twobodykeyword
        self.kbodyfile = kbodyfile
        self._fold = None

    def run(self):
        self._readbond()
        self._rungaussian()
        self._takeforce()

    def _rungaussian(self):
        if not self.gaussiancommand:
            GaussianRunner(command=self.command, cpu_num=self.nproc_sum, nproc=self.nproc).runGaussianInParallel(
                'gjf', [os.path.join(self.gaussian_dir, job+".gjf") for job in self.jobs])
        else:
            with open(self.jobfile, 'w') as f:
                print(*[os.path.join(self.gaussian_dir, job+".gjf")
                        for job in self.jobs], file=f)
            sp.call(self.gaussiancommand.split())

    def _logging(self, *message):
        localtime = time.asctime(time.localtime(time.time()))
        print(localtime, 'AIMDFragmentation', *message)

    def _getjobname(self, *molid):
        if len(molid) == 1:
            return f"mol{molid[0]}"
        elif len(molid) == 2:
            molid = sorted(molid)
            return f'tb{molid[0]}-{molid[1]}'

    def _mo(self, i, bond, molecule, done):  # connect molecule
        molecule.append(i)
        done[i] = True
        for b in bond[i]:
            if not done[b]:
                molecule, done = self._mo(b, bond, molecule, done)
        return molecule, done

    def _readpdb(self):
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats('xyz', 'pdb')
        mol = openbabel.OBMol()
        conv.ReadFile(mol, self.xyzfilename)
        pdbstring = conv.WriteString(mol)
        bond = [[] for x in range(self._natom)]
        for line in pdbstring.split('\n'):
            if line.startswith("CONECT"):
                s = line.split()
                bond[int(s[1])-1] += [int(x)-1 for x in s[2:]]
        d, done = [], [False for x in range(self._natom)]
        for i in range(self._natom):
            if not done[i]:
                molecule, done = self._mo(i, bond, [], done)
                d.append(molecule)
        self._mols = d

    def _printgjf(self, jobname, selected_atomsid):
        if not os.path.exists(self.gaussian_dir):
            os.makedirs(self.gaussian_dir)
        buff = []
        # only supported for C, H, and O
        S = ((3 if atoms.get_chemical_symbols() == ["O", "O"] else (
            Counter(atoms.get_chemical_symbols())['H'] % 2+1)) for atoms in [self._atoms[atomsid] for atomsid in selected_atomsid])
        atoms_whole, S_whole, kbodykeyword = self._atoms[sum(selected_atomsid, [])], sum(
            S)-len(selected_atomsid)+1, (self.onebodykeyword if len(selected_atomsid) == 1 else self.twobodykeyword)
        atoms_whole.wrap(
            center=atoms_whole[0].position/atoms_whole.get_cell_lengths_and_angles()[0:3], pbc=atoms_whole.get_pbc())
        buff.append(
            f'%nproc={self.nproc}\n%mem={self.qmmem}\n# force {self.qmmethod}/{self.qmbasis} {kbodykeyword} {self.addkw}\n\n{jobname}\n\n0 {S_whole}')
        buff.extend(("{} {} {} {}".format(atom.symbol, *atom.position)
                     for atom in atoms_whole))
        buff.append('\n')
        with open(os.path.join(self.gaussian_dir, jobname+".gjf"), 'w') as f:
            f.write('\n'.join(buff))

    def _processjob(self, molid, selected_atoms):
        jobname = self._getjobname(*molid)
        self._atomid[jobname] = sum(selected_atoms, [])
        self._printgjf(jobname, selected_atoms)
        self.jobs.append(jobname)

    def _printmol(self):
        for molid, atoms in enumerate(self._mols, 1):
            self._processjob((molid,), (atoms,))

    def _printtb(self):
        for molid1, atoms1 in enumerate(self._mols, 1):
            for molid2, atoms2 in enumerate(self._mols[molid1:], molid1+1):
                if np.min(get_distances(self._atoms[atoms1].positions, self._atoms[atoms2].positions, cell=self._atoms.get_cell(), pbc=self._atoms.get_pbc())[1]) <= self.cutoff:
                    self._processjob((molid1, molid2), (atoms1, atoms2))

    def _readbond(self):
        self._atoms = readxyz(self.xyzfilename)
        self._atoms.set_pbc(self.pbc)
        self._atoms.set_cell(self.cell)
        self._natom = len(self._atoms)
        self._readpdb()
        self._printmol()
        self._printtb()

    def _readforce(self, jobname):
        forces = GaussianAnalyst(properties=['force']).readFromLOG(
            os.path.join(self.gaussian_dir, jobname+'.log'))['force']
        if forces:
            atoms = {}
            for index, force in forces.items():
                atoms[self._atomid[jobname][index-1]
                      ] = np.array(force)*self.unit
            return atoms, None
        else:
            return None, jobname

    @property
    def fold(self):
        if self._fold is None:
            if os.path.isfile(self.kbodyfile):
                loadingfold = np.loadtxt(self.kbodyfile)
                if loadingfold.shape == (self._natom, 6):
                    self._fold = loadingfold
                    self._logging("Load old forces.")
        if self._fold is None:
            self._fold = np.zeros((self._natom, 6))
            self._logging("No old forces found. Use 0 instead.")
        return self._fold

    def _takeforce(self):
        onebodyforce, twobodyforce = np.zeros(
            (self._natom, 3)), np.zeros((self._natom, 3))
        with Pool(self.nproc_sum) as pool:
            onebodyresults = pool.imap(
                self._readforce, [self._getjobname(i) for i in range(1, len(self._mols)+1)])
            twobodyresults = pool.imap(self._readforce, [self._getjobname(i, j) for i in range(1, len(
                self._mols)+1) for j in range(i+1, len(self._mols)+1) if self._getjobname(i, j) in self.jobs])
            twobodyerroratoms = []
            for i, results in enumerate((onebodyresults, twobodyresults)):
                for atoms, jobname in results:
                    if atoms:
                        for atom, force in atoms.items():
                            if i == 0:
                                onebodyforce[atom] = force
                            else:
                                twobodyforce[atom] += force-onebodyforce[atom]
                    else:
                        self._logging('WARNING:', 'No forces of', jobname,
                                      'found. Use the old forces instead.')
                        self.errorfiles.append(os.path.join(
                            self.gaussian_dir, jobname+'.log'))
                        if i == 0:
                            onebodyforce[self._atomid[jobname]
                                         ] = self.fold[self._atomid[jobname]][:, 0:3]
                        else:
                            twobodyerroratoms += self._atomid[jobname]
            if twobodyerroratoms:
                twobodyforce[twobodyerroratoms] = self.fold[twobodyerroratoms][:, 3:6]
                self._logging("Atom", *twobodyerroratoms,
                              "use(s) the old 2-body forces.")
        finalforces = onebodyforce+twobodyforce
        # Make the resultant force equal to 0
        if np.abs(np.sum(finalforces)) > 0:
            finalforces -= np.abs(finalforces) / \
                np.sum(np.abs(finalforces), 0)*np.sum(finalforces, 0)
        np.savetxt(self.outputfile, finalforces, fmt='%16.9f')
        np.savetxt(self.kbodyfile, np.hstack(
            (onebodyforce, twobodyforce)), fmt='%16.9f')
        forcesum = np.sum(finalforces, axis=0)
        forcesumdis = np.linalg.norm(forcesum)
        self._logging("Resultant force:", *("%16.9f" % x for x in forcesum))
        self._logging("Magnitude:", "%16.9f" % forcesumdis)
