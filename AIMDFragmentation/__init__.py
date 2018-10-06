#!/usr/bin/env python3
# updated at 2018/7/23 3:00
# Add root of readbond.py to PATH
# Usage:
# readbond.py 28 xxx.xyz
import numpy as np
import sys
import os
from GaussianRunner import GaussianRunner,GaussianAnalyst

class AIMDFragmentation(object):
    def __init__(self,nproc_sum,nproc,cutoff,xyzfilename,pdbfilename,qmmethod,qmbasis,addkw,qmmem,atombondnumber,logfile,outputfile="force.dat",unit=1):
        self.nproc_sum=nproc_sum
        self.nproc=nproc
        self.cutoff=cutoff
        self.xyzfilename=xyzfilename
        self.pdbfilename=pdbfilename
        self.qmmethod=qmmethod
        self.qmbasis=qmbasis
        self.addkw=addkw
        self.qmmem=qmmem
        self.atombondnumber=atombondnumber
        self.logfile=logfile
        self.outputfile=outputfile
        self.unit=unit

    def run(self):
        os.system('obabel -ixyz '+self.xyzfilename+' -opdb -O '+self.pdbfilename+' > /dev/null')
        g16commands=self.readbond()
        GaussianRunner(command='g16',cpu_num=self.nproc_sum,nproc=self.nproc).runGaussianInParallel('gjf',g16commands)
        self.takeforce(g16=g16commands)

    def mo(self,i,bond,molecule,done,bondlist): #connect molecule
        molecule.append(i)
        done[i]=True
        for j in range(len(bond[i])):
            b=bond[i][j]
            bo=(i,b,1) if i<b else (b,i,1)
            if not bo in bondlist:
                bondlist.append(bo)
            if not b in done:
                molecule,done,bondlist=self.mo(b,bond,molecule,done,bondlist)
        return molecule,done,bondlist

    def printmoleculename(self,atoms,atomtype):
        typenumber={}
        for atomnumber in atoms:
            if atomtype[atomnumber] in typenumber:
                typenumber[atomtype[atomnumber]]+=1
            else:
                typenumber[atomtype[atomnumber]]=1
        name="".join([key+(str(value) if value>1 else "") for key,value in typenumber.items()])
        return name

    def readpdb(self):
        bond={}
        atomtype={}
        atomxyz={}
        with open(self.pdbfilename) as f:
            for line in f:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    s=line.split()
                    atomtype[int(s[1])]=s[2]
                    atomxyz[int(s[1])]=np.array([float(s[5]),float(s[6]),float(s[7])])
                    bond[int(s[1])]=[]
                elif line.startswith("CONECT"):
                    s=line.split()
                    bond[int(s[1])]+=[int(x) for x in s[2:]]
        d={}
        done={}
        molid=0
        for i in bond.keys():
            if not i in done:
                mole,done,bondlist=self.mo(i,bond,[],done,[])
                mole.sort()
                bondlist.sort()
                molid+=1
                d[molid]=(self.printmoleculename(mole,atomtype),mole,bondlist)
        return d,atomtype,atomxyz

    def printgjf(self,jobname,xyzoutput,S):
        with open(jobname+".gjf",'w') as f:
            print("%nproc="+str(self.nproc)+"\n%mem="+self.qmmem+"\n# force "+self.qmmethod+"/"+self.qmbasis+" "+self.addkw+"\n\n"+jobname+"\n\n0 "+str(S),file=f)
            print("\n".join(xyzoutput),file=f)
            print("",file=f)

    def getxyz(self,atoms,atomtype,atomxyz,jobname):
        atomindex={}
        xyzoutput=[]
        atomnumber=0
        with open(jobname+".id",'w') as f:
            for atom in atoms:
                atomnumber+=1
                atomindex[atom]=atomnumber
                print(atomnumber,atom,file=f)
                type=atomtype[atom]
                x,y,z=atomxyz[atom]
                xyzoutput.append(" "+type+" "+str(x)+" "+str(y)+" "+str(z))
        return xyzoutput,atomindex,atomnumber

    def printmol(self,d,atomtype,atomxyz):
        Smol={}
        g16=[]
        for molid,(moleculename,atoms,bonds) in d.items():
            jobname="mol"+str(molid)
            xyzoutput,atomindex,atomnumber=self.getxyz(atoms,atomtype,atomxyz,jobname)
            bondnumber=[0 for x in range(atomnumber)]
            for bond in bonds:
                for i in range(2):
                    bondnumber[atomindex[bond[i]]-1]+=bond[2]
            S=1
            for atom,index in atomindex.items():
                S+=self.atombondnumber[atomtype[atom]]-bondnumber[index-1]
            S= 3 if moleculename=="O2" else (2 if S%2==0 else 1)
            Smol[molid]=S
            self.printgjf(jobname,xyzoutput,S)
            g16.append(jobname+".gjf")
        return Smol,g16

    def printtb(self,d,atomtype,atomxyz,Smol):
        g16=[]
        for molid1,(moleculename1,atoms1,bonds1) in d.items():
            for molid2,(moleculename2,atoms2,bonds2) in d.items():
                if molid1>=molid2:
                    continue
                for atom1 in atoms1:
                    for atom2 in atoms2:
                        if np.linalg.norm(atomxyz[atom1]-atomxyz[atom2])<=self.cutoff:
                            break
                    else:
                        continue
                    break
                else:
                    continue
                atoms=atoms1+atoms2
                bonds=bonds1+bonds2
                jobname="tb"+str(molid1)+"-"+str(molid2)
                xyzoutput,atomindex,atomnumber=self.getxyz(atoms,atomtype,atomxyz,jobname)
                S=Smol[molid1]+Smol[molid2]-1
                S=2 if S%2==0 else 1
                for molid in (molid1,molid2):
                    if Smol[molid]==3:
                        S+=2
                self.printgjf(jobname,xyzoutput,S)
                g16.append(jobname+".gjf")
        return g16

    def readbond(self):
        d,atomtype,atomxyz=self.readpdb()
        Smol,g16mol=self.printmol(d,atomtype,atomxyz)
        with open(self.logfile,'a') as f:
            print("Total S",sum(Smol.values())-len(Smol)+1,file=f)
        g16tb=self.printtb(d,atomtype,atomxyz,Smol)
        return g16mol+g16tb

    def readforce(self,jobname):
        d={}
        with open(jobname+".id") as f:
            for line in f:
                index,atomid=(int(x) for x in line.split())
                d[index]=atomid
        forces=GaussianAnalyst(properties=['force']).readFromLOG(jobname+'.log')['force']
        atoms={}
        for index,force in forces.items():
            atoms[d[index]]=np.array(force)
        return atoms

    def takeforce(self,g16):
        atomforce={}
        i=0
        while "mol"+str(i+1)+".gjf" in g16:
            i+=1
            forces=self.readforce("mol"+str(i))
            for atom,force in forces.items():
                atomforce[atom]=force
        n=i
        twobodyforce={}
        for i in range(1,n+1):
            for j in range(i+1,n+1):
                if not ("tb"+str(i)+"-"+str(j)+".gjf" in g16):
                    continue
                forces=self.readforce("tb"+str(i)+"-"+str(j))
                for atom,force in forces.items():
                    twobodyforce[atom]=twobodyforce[atom]+force-atomforce[atom] if atom in twobodyforce else force-atomforce[atom]
        finalforces=[]
        with open(self.outputfile,'w') as f:
            for atom,force in sorted(atomforce.items(),key=lambda item:item[0]):
                finalforce=force+twobodyforce[atom] if atom in twobodyforce else force
                finalforce*=self.unit
                print ("".join("%16.9f"%x for x in finalforce),file=f)
                finalforces.append(finalforce)
        with open(self.logfile,'a') as f:
            forcesum=np.sum(finalforces,axis=0)
            print("".join("%16.9f"%x for x in forcesum),file=f,end='')
            forcesumdis=np.sqrt(np.sum(np.square(forcesum)))
            print("%16.9f"%forcesumdis,file=f)
            print("\n",file=f)

AIMDBlock=AIMDFragmentation

def readmfccin(mfccinfilename="mfcc.in"):
    glb={}
    loc={}
    with open(mfccinfilename) as f:
        for line in f:
            line=line.strip('\n')
            if line=="&qm_ctrl":
                qmctrl=True
            elif qmctrl:
                if line=="&end":
                    break
                else:
                    exec(line.strip(),glb,loc)
    qmproc,qmmethod,qmbasis,addkw,dl,qmmem=loc['qmproc'],loc['qmmethod'],loc['qmbasis'],loc['addkw'],loc['dl'],loc['qmmem']
    return qmproc,qmmethod,qmbasis,addkw,dl,qmmem

if __name__ == '__main__':
    xyzfilename=sys.argv[2]
    pdbfilename=xyzfilename+".pdb"
    atombondnumber={"C":4,"H":1,"O":2}
    qmproc,qmmethod,qmbasis,addkw,dl,qmmem=readmfccin(mfccinfilename="mfcc.in")
    AIMDFragmentation(nproc_sum=int(sys.argv[1]),nproc=qmproc,cutoff=dl,xyzfilename=xyzfilename,pdbfilename=pdbfilename,qmmethod=qmmethod,qmbasis=qmbasis,addkw=addkw,qmmem=qmmem,atombondnumber=atombondnumber,logfile="force.log").run()
