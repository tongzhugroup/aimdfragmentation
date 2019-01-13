#!/usr/bin/env python3
"""Fragmentation for AIMD"""

__author__ = "Jinzhe Zeng"
__email__ = "jzzeng@stu.ecnu.edu.cn"
__update__ = "2019-01-13"
__date__ = "2018-07-18"


import itertools
import os
import subprocess as sp
import sys
import time
from collections import Counter, defaultdict
from functools import partial
from multiprocessing import Pool, cpu_count

import numpy as np
import openbabel
from ase import Atoms
from ase.geometry import get_distances
from ase.io import read as readxyz
from gaussianrunner import GaussianAnalyst, GaussianRunner


class AIMDFragmentation(object):
    def __init__(self, nproc_sum=None, nproc=4, cutoff=3.5, xyzfilename="comb.xyz", pdbfilename="comb.pdb", qmmethod="mn15", qmbasis="6-31g(d)", addkw="", qmmem="400MW", atombondnumber=None, logfile=None, outputfile="force.dat", outputenergyfile="energy.dat", unit=1, energyunit=1, pbc=False, cell=[0, 0, 0], gaussian_dir="gaussian_files", command="g16", gaussiancommand=None, jobfile="gaussianjobs", onebodykeyword="scf=(xqc,MaxConventionalCycles=256)", twobodykeyword="scf=(maxcyc=256)", kbodyfile="kbforce.dat", fg=True, kmax=3):
        self.nproc_sum = nproc_sum if nproc_sum else cpu_count()
        self.nproc = nproc
        self.cutoff = cutoff
        self.xyzfilename = xyzfilename
        self.qmmethod = qmmethod
        self.qmbasis = qmbasis
        self.addkw = addkw
        self.qmmem = qmmem
        self.outputfile = outputfile
        self.outputenergyfile = outputenergyfile
        self.unit = unit
        self.energyunit = energyunit
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
        self.fg = fg
        self.kmax = kmax
        self._distances = {}

    def run(self):
        self._readbond()
        self._rungaussian()
        self._takeforce()

    def _rungaussian(self):
        if not self.gaussiancommand:
            GaussianRunner(command=self.command, cpu_num=self.nproc_sum, nproc=self.nproc).runGaussianInParallel(
                'gjf', [os.path.join(self.gaussian_dir, f"{job}.gjf") for job in self.jobs])
        else:
            with open(self.jobfile, 'w') as f:
                print(*[os.path.join(self.gaussian_dir, f"{job}.gjf")
                        for job in self.jobs], file=f)
            sp.call(self.gaussiancommand.split())

    def _logging(self, *message):
        localtime = time.asctime(time.localtime(time.time()))
        print(localtime, 'AIMDFragmentation', *message)

    def _getjobname(self, *molid):
        molid = sorted(molid)
        return f'{len(molid)}b{"-".join((str(x) for x in molid))}'

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
        selected_atoms = [self._atoms[atomsid] for atomsid in selected_atomsid]
        multiplicities = list((3 if atoms.get_chemical_symbols() == ["O", "O"] else (
            Counter(atoms.get_chemical_symbols())['H'] % 2+1)) for atoms in selected_atoms)
        atoms_whole, multiplicity_whole, kbodykeyword = self._atoms[sum(selected_atomsid, [])], sum(
            multiplicities)-len(selected_atomsid)+1, (self.onebodykeyword if len(selected_atomsid) == 1 else self.twobodykeyword)
        nproc = f'%nproc={self.nproc}'
        title = f'\n{jobname}\n'
        mem = f'mem={self.qmmem}'
        if len(selected_atomsid) == 1 or not self.fg:
            atoms_whole.wrap(
                center=atoms_whole[0].position/atoms_whole.get_cell_lengths_and_angles()[0:3], pbc=atoms_whole.get_pbc())
            kw = f'# force {self.qmmethod}/{self.qmbasis} {kbodykeyword} {self.addkw}'
            buff.extend((nproc, mem, kw, title))
            buff.append(f'0 {multiplicity_whole}')
            buff.extend(("{} {} {} {}".format(atom.symbol, *atom.position)
                         for atom in atoms_whole))
            buff.append('\n')
        else:
            for atoms in selected_atoms:
                atoms.wrap(center=atoms_whole[0].position/atoms_whole.get_cell_lengths_and_angles()[
                           0:3], pbc=atoms_whole.get_pbc())
            chk = f'%chk={os.path.join(self.gaussian_dir, jobname+".chk")}'
            connect = '\n--link1--\n'
            kw1 = f'# {self.qmmethod}/{self.qmbasis} {kbodykeyword} {self.addkw} guess=fragment={len(selected_atomsid)}'
            kw2 = f'# force {self.qmmethod}/{self.qmbasis} {kbodykeyword} {self.addkw} geom=chk guess=read'
            multiplicities_str = ' '.join(
                (f'0 {multiplicity}' for multiplicity in itertools.chain((multiplicity_whole,), multiplicities)))
            buff.extend((chk, nproc, mem, kw1, title, multiplicities_str))
            for index, atoms in enumerate(selected_atoms, 1):
                buff.extend(('{}(Fragment={}) {} {} {}'.format(
                    atom.symbol, index, *atom.position) for atom in atoms))
            buff.extend((connect, chk, nproc, kw2,
                         title, multiplicities_str, '\n'))
        with open(os.path.join(self.gaussian_dir, f"{jobname}.gjf"), 'w') as f:
            f.write('\n'.join(buff))

    def _processjob(self, molids):
        jobname = self._getjobname(*molids)
        selected_atoms = [self._mols[molid-1] for molid in molids]
        self._atomid[jobname] = sum(selected_atoms, [])
        self._printgjf(jobname, selected_atoms)
        self.jobs.append(jobname)

    def _printkb(self, k):
        for molids in itertools.combinations(range(1, len(self._mols)+1), k):
            if all((self._isclose(molida, molidb) for molida, molidb in itertools.combinations(molids, 2))):
                self._processjob(molids)

    def _isclose(self, molid1, molid2):
        name = f'{molid1}-{molid2}'
        if not name in self._distances:
            self._distances[name] = np.min(get_distances(self._atoms[self._mols[molid1-1]].positions, self._atoms[self._mols[molid2-1]].positions,
                                                         cell=self._atoms.get_cell(), pbc=self._atoms.get_pbc())[1]) <= self.cutoff
        return self._distances[name]

    def _readbond(self):
        self._atoms = readxyz(self.xyzfilename)
        self._atoms.set_pbc(self.pbc)
        self._atoms.set_cell(self.cell)
        self._natom = len(self._atoms)
        self._readpdb()
        for k in range(1, self.kmax+1):
            self._printkb(k)

    def _readforce(self, mols):
        jobname = self._getjobname(*mols)
        results = GaussianAnalyst(properties=['force', 'energy']).readFromLOG(
            os.path.join(self.gaussian_dir, f'{jobname}.log'))
        qm_force = results['force']
        qm_energy = results['energy']
        if not qm_force is None:
            forces = np.zeros((self._natom, 3))
            forces[self._atomid[jobname]] = qm_force
            forces *= self.unit
            energy = qm_energy * self.energyunit
            return forces, energy, mols
        else:
            return None, None, mols

    @property
    def fold(self):
        if self._fold is None:
            if os.path.isfile(self.kbodyfile):
                loadingfold = np.loadtxt(self.kbodyfile)
                if loadingfold.shape == (self._natom, 3*self.kmax):
                    self._fold = loadingfold
                    self._logging("Load old forces.")
        if self._fold is None:
            self._fold = np.zeros((self._natom, 3*self.kmax))
            self._logging("No old forces found. Use 0 instead.")
        return self._fold

    def _takeforce(self):
        kbodyforces = np.zeros((self.kmax, self._natom, 3))
        kbodyenergies = np.zeros((self.kmax,))
        kbodyerroratoms = [list() for x in range(self.kmax)]
        molsforces = defaultdict(partial(np.zeros, (self._natom, 3)))
        molsenergies = defaultdict(float)
        with Pool(self.nproc_sum) as pool:
            kbodyresults = [pool.imap(
                self._readforce, [molids for molids in itertools.combinations(range(1, len(self._mols)+1), i+1) if self._getjobname(*molids) in self.jobs]) for i in range(self.kmax)]
            for i, results in enumerate(kbodyresults):
                for force, energy, mols in results:
                    if not force is None:
                        molsforces[mols] = force - np.sum((np.sum(
                            (molsforces[klessmols] for klessmols in itertools.combinations(mols, j)), axis=0) for j in range(i+1)), axis=0)
                        molsenergies[mols] = energy - np.sum((np.sum(
                            (molsenergies[klessmols] for klessmols in itertools.combinations(mols, j))) for j in range(i+1)))
                        kbodyforces[i] += molsforces[mols]
                        kbodyenergies[i] += molsenergies[mols]
                    else:
                        jobname = self._getjobname(*mols)
                        self._logging(
                            f'WARNING: No forces of {jobname} found. Use the old forces instead.')
                        self.errorfiles.append(os.path.join(
                            self.gaussian_dir, f'{jobname}.log'))
                        if i == 0:
                            kbodyforces[i][self._atomid[jobname]
                                           ] = self.fold[self._atomid[jobname]][:, 0:3]
                        else:
                            kbodyerroratoms[i] += self._atomid[jobname]
            for i in range(1, self.kmax):
                if kbodyerroratoms[i]:
                    kbodyforces[i][kbodyerroratoms[i]
                                   ] = self.fold[kbodyerroratoms[i]][:, 3:6]
                    self._logging("Atom", *kbodyerroratoms[i],
                                  f"use(s) the old {i}-body forces.")
        finalforces = np.sum(kbodyforces, axis=0)
        finalenergies = np.sum(kbodyenergies)
        # Make the resultant force equal to 0
        if np.abs(np.sum(finalforces)) > 0:
            finalforces -= np.abs(finalforces) / \
                np.sum(np.abs(finalforces), 0)*np.sum(finalforces, 0)
        np.savetxt(self.outputfile, finalforces, fmt='%16.9f')
        with open(self.outputenergyfile, 'w') as f:
            f.write(f'{finalenergies:16.9f}')
        np.savetxt(self.kbodyfile, np.hstack(
            (kbodyforces[i] for i in range(self.kmax))), fmt='%16.9f')
        forcesum = np.sum(finalforces, axis=0)
        forcesumdis = np.linalg.norm(forcesum)
        self._logging(f"Energy: {finalenergies:16.9f}")
        self._logging("Resultant force:", *(f"{x:16.9f}" for x in forcesum))
        self._logging(f"Magnitude: {forcesumdis:16.9f}")
