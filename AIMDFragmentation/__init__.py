#!/usr/bin/env python3
import numpy as np
import sys, os, time
from GaussianRunner import GaussianRunner,GaussianAnalyst
from ase.io import read as readxyz
from ase.geometry import get_distances

class AIMDFragmentation(object):
    def __init__(self,nproc_sum,nproc,cutoff,xyzfilename,pdbfilename,qmmethod,qmbasis,addkw,qmmem,atombondnumber={"C":4,"H":1,"O":2},logfile="force.log",outputfile="force.dat",unit=1,pbc=False,cell=[0,0,0],gaussian_dir="gaussian_files"):
        self.nproc_sum=nproc_sum
        self.nproc=nproc
        self.cutoff=cutoff
        self.xyzfilename=xyzfilename
        self.pdbfilename=pdbfilename
        self.qmmethod=qmmethod
        self.qmbasis=qmbasis
        self.addkw=addkw
        self.qmmem=qmmem
        #self.atombondnumber=atombondnumber
        self.logfile=logfile
        self.outputfile=outputfile
        self.unit=unit
        self.openlogfile=False
        self.pbc=pbc
        self.cell=cell
        self.atomid={}
        self.jobs=[]
        self.gaussian_dir=gaussian_dir

    def run(self):
        self.readbond()
        GaussianRunner(command='g16',cpu_num=self.nproc_sum,nproc=self.nproc).runGaussianInParallel('gjf',[os.path.join(self.gaussian_dir,job+".gjf") for job in self.jobs])
        self.takeforce()

    def logging(self,*message):
        if not self.openlogfile:
            self.openlogfile=open(self.logfile,'a')
        localtime = time.asctime( time.localtime(time.time()) )
        print(localtime,'AIMDFragmentation',*message)
        print(localtime,'AIMDFragmentation',*message,file=self.openlogfile)

    def getjobname(self,*molid):
        if len(molid)==1:
            return "mol"+str(*molid)
        elif len(molid)==2:
            molid1,molid2=sorted(molid)
            return "tb"+str(molid1)+"-"+str(molid2)

    def mo(self,i,bond,molecule,done): #connect molecule
        molecule.append(i)
        done[i]=True
        for b in bond[i]:
            if not done[b]:
                molecule,done=self.mo(b,bond,molecule,done)
        return molecule,done

    def readpdb(self):
        bond=[[] for x in range(self.natom)]
        with open(self.pdbfilename) as f:
            for line in f:
                if line.startswith("CONECT"):
                    s=line.split()
                    bond[int(s[1])-1]+=[int(x)-1 for x in s[2:]]
        d=[]
        done=[False for x in range(self.natom)]
        for i in range(self.natom):
            if not done[i]:
                molecule,done=self.mo(i,bond,[],done)
                molecule.sort()
                d.append(molecule)
        self.mols=d

    def printgjf(self,jobname,selected_atoms1,S1,selected_atoms2=None,S2=None):
        selected_atoms1.wrap(center=selected_atoms1[0].position/selected_atoms1.get_cell_lengths_and_angles()[0:3],pbc=selected_atoms1.get_pbc())
        if not os.path.exists(self.gaussian_dir):
            os.makedirs(self.gaussian_dir)
        with open(os.path.join(self.gaussian_dir,jobname+".gjf"),'w') as f:
            if not selected_atoms2:
                print("%nproc="+str(self.nproc),file=f)
                print("%mem="+self.qmmem,file=f)
                print("#","force",self.qmmethod+"/"+self.qmbasis,self.addkw+"\n\n"+jobname+"\n\n0",S1,file=f)
                for atom in selected_atoms1:
                    print(atom.symbol,*atom.position,file=f)
                print("",file=f)
            else:
                if S1>S2:
                    S2*=-1
                else:
                    S1*=-1
                Stotal=S1+S2+1
                print("%chk="+os.path.join(self.gaussian_dir,jobname+".chk"),file=f)
                print("%nproc="+str(self.nproc),file=f)
                print("%mem="+self.qmmem,file=f)
                print("#",self.qmmethod+"/"+self.qmbasis,"guess=fragment=2",self.addkw+"\n\n"+jobname+"\n\n0",Stotal,0,S1,0,S2,file=f)
                for index,selected_atoms in enumerate((selected_atoms1,selected_atoms2),start=1):
                    for atom in selected_atoms:
                        print(atom.symbol+"(Fragment="+str(index)+")",*atom.position,file=f)
                print("",file=f)
                print("--link1--",file=f)
                print("%chk="+os.path.join(self.gaussian_dir,jobname+".chk"),file=f)
                print("%nproc="+str(self.nproc),file=f)
                print("%mem="+self.qmmem,file=f)
                print("#",self.qmmethod+"/"+self.qmbasis,"guess=read","geom=chk","force",self.addkw+"\n\n"+jobname+"\n\n0",Stotal,0,S1,0,S2,file=f)
                print("",file=f)

    def printmol(self):
        self.Smol=[]
        for molid,atoms in enumerate(self.mols,1):
            jobname=self.getjobname(molid)
            self.atomid[jobname]=atoms
            selected_atoms=self.atoms[atoms]
            # only supported for C,H,O
            num_H=sum(1 for atom in selected_atoms if atom.symbol=="H")
            S= 3 if selected_atoms.get_chemical_symbols()==['O','O'] else (1 + num_H % 2)
            self.Smol.append(S)
            self.printgjf(jobname,selected_atoms,S)
            self.jobs.append(jobname)

    def printtb(self):
        for molid1,atoms1 in enumerate(self.mols,1):
            for molid2,atoms2 in enumerate(self.mols[molid1:],molid1+1):
                if np.min(get_distances(self.atoms[atoms1].positions,self.atoms[atoms2].positions,cell=self.atoms.get_cell(),pbc=self.atoms.get_pbc())[1])<=self.cutoff:
                    jobname=self.getjobname(molid1,molid2)
                    self.atomid[jobname]=atoms1+atoms2
                    selected_atoms1,selected_atoms2=self.atoms[atoms1],self.atoms[atoms2]
                    S1,S2=self.Smol[molid1-1],self.Smol[molid2-1]
                    self.printgjf(jobname,selected_atoms1,S1,selected_atoms2,S2)
                    self.jobs.append(jobname)

    def readbond(self):
        self.atoms=readxyz(self.xyzfilename)
        self.atoms.set_pbc(self.pbc)
        self.atoms.set_cell(self.cell)
        self.natom=len(self.atoms)
        os.system('obabel -ixyz '+self.xyzfilename+' -opdb -O '+self.pdbfilename+' > /dev/null')
        self.readpdb()
        self.printmol()
        self.logging("Total S",sum(self.Smol)-len(self.Smol)+1)
        self.printtb()

    def readforce(self,jobname):
        forces=GaussianAnalyst(properties=['force']).readFromLOG(os.path.join(self.gaussian_dir,jobname+'.log'))['force']
        atoms={}
        for index,force in forces.items():
            atoms[self.atomid[jobname][index-1]]=np.array(force)
        return atoms

    def takeforce(self):
        onebodyforce=np.zeros((self.natom,3))
        for i in range(1,len(self.mols)+1):
            forces=self.readforce(self.getjobname(i))
            for atom,force in forces.items():
                onebodyforce[atom]=force
        twobodyforce=np.zeros((self.natom,3))
        for i in range(1,len(self.mols)+1):
            for j in range(i+1,len(self.mols)+1):
                if self.getjobname(i,j) in self.jobs:
                    forces=self.readforce("tb"+str(i)+"-"+str(j))
                    for atom,force in forces.items():
                        twobodyforce[atom]=twobodyforce[atom]+force-onebodyforce[atom]
        finalforces=onebodyforce+twobodyforce
        finalforces*=self.unit
        with open(self.outputfile,'w') as f:
            for finalforce in finalforces:
                print(*("%16.9f"%x for x in finalforce),sep="",file=f)
        forcesum=np.sum(finalforces,axis=0)
        self.logging("Total force:",*("%16.9f"%x for x in forcesum))
        forcesumdis=np.linalg.norm(forcesum)
        self.logging("Magnitude:","%16.9f"%forcesumdis)
