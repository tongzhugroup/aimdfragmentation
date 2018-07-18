#!/usr/bin/env python3
# updated at 2018/5/31 17:00
# Add root of readbond.py to PATH
# Usage:
# readbond.py 28 xxx.xyz
import numpy as np
import sys
import os
def mo(i,bond,molecule,done,bondlist): #connect molecule
    molecule.append(i)
    done[i]=True
    for j in range(len(bond[i])):
        b=bond[i][j]
        bo=(i,b,1) if i<b else (b,i,1)
        if not bo in bondlist:
            bondlist.append(bo)
        if not b in done:
            molecule,done,bondlist=mo(b,bond,molecule,done,bondlist)
    return molecule,done,bondlist

def printmoleculename(atoms,atomtype):
    typenumber={}
    for atomnumber in atoms:
        if atomtype[atomnumber] in typenumber:
            typenumber[atomtype[atomnumber]]+=1
        else:
            typenumber[atomtype[atomnumber]]=1
    name="".join([key+(str(value) if value>1 else "") for key,value in typenumber.items()])                
    return name
    
def readpdb(pdbfilename):
    bond={}
    atomtype={}
    atomxyz={}
    with open(pdbfilename) as f:
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
            mole,done,bondlist=mo(i,bond,[],done,[])
            mole.sort()
            bondlist.sort()
            molid+=1
            d[molid]=(printmoleculename(mole,atomtype),mole,bondlist)
    return d,atomtype,atomxyz

def printgjf(nproc,jobname,xyzoutput,bondoutput,S,qmmethod,qmbasis,addkw,qmmem):
    with open(jobname+".gjf",'w') as f:
        print("%nproc="+str(nproc)+"\n%mem="+qmmem+"\n# force "+qmmethod+"/"+qmbasis+" "+addkw+" geom=connectivity\n\n"+jobname+"\n\n0 "+str(S),file=f)
        print("\n".join(xyzoutput),file=f)
        print("",file=f)
        print("\n".join(bondoutput),file=f)
        print("",file=f)
    return 
    
def getxyz(atoms,atomtype,atomxyz,jobname):
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
def printmol(nproc,d,atomtype,atomxyz,atombondnumber,qmmethod,qmbasis,addkw,qmmem):
    Smol={}
    g16=[]
    for molid,molecule in d.items():
        moleculename,atoms,bonds=molecule
        jobname="mol"+str(molid)
        xyzoutput,atomindex,atomnumber=getxyz(atoms,atomtype,atomxyz,jobname)
        bondoutput=[" "+str(x+1) for x in range(atomnumber)]
        bondnumber=[0 for x in range(atomnumber)]
        for bond in bonds:
            bondoutput[atomindex[bond[0]]-1]+=" "+str(atomindex[bond[1]])+" "+str(bond[2])
            bondnumber[atomindex[bond[0]]-1]+=bond[2]
            bondnumber[atomindex[bond[1]]-1]+=bond[2]
        S=1
        for atom,index in atomindex.items():
            S+=atombondnumber[atomtype[atom]]-bondnumber[index-1]
        S= 3 if moleculename=="O2" else (2 if S%2==0 else 1)
        Smol[molid]=S
        printgjf(nproc,jobname,xyzoutput,bondoutput,S,qmmethod,qmbasis,addkw,qmmem)
        g16.append("g16 "+jobname+".gjf")
    return Smol,g16

   
def printtb(nproc,d,atomtype,atomxyz,Smol,r,qmmethod,qmbasis,addkw,qmmem):
    g16=[]
    for molid1,molecule1 in d.items():
        moleculename1,atoms1,bonds1=molecule1
        for molid2,molecule2 in d.items():
            moleculename2,atoms2,bonds2=molecule2
            if molid1>=molid2:
                continue
            for atom1 in atoms1:
                for atom2 in atoms2:
                    if np.sum(np.square(atomxyz[atom1]-atomxyz[atom2]))<=r**2:
                        break
                else:
                    continue
                break
            else:
                continue
            atoms=atoms1+atoms2
            bonds=bonds1+bonds2
            jobname="tb"+str(molid1)+"-"+str(molid2)
            xyzoutput,atomindex,atomnumber=getxyz(atoms,atomtype,atomxyz,jobname)
            bondoutput=[" "+str(x+1) for x in range(atomnumber)]
            for bond in bonds:
                bondoutput[atomindex[bond[0]]-1]+=" "+str(atomindex[bond[1]])+" "+str(bond[2])
            S=Smol[molid1]+Smol[molid2]-1
            S=2 if S%2==0 else 1
            for molid in (molid1,molid2):
                if Smol[molid]==3:
                    S+=2
            printgjf(nproc,jobname,xyzoutput,bondoutput,S,qmmethod,qmbasis,addkw,qmmem)
            g16.append("g16 "+jobname+".gjf")
    return g16
            
def readbond(nproc,pdbfilename,r,qmmethod,qmbasis,addkw,qmmem,atombondnumber,logfile):
    d,atomtype,atomxyz=readpdb(pdbfilename)
    Smol,g16mol=printmol(nproc,d,atomtype,atomxyz,atombondnumber,qmmethod=qmmethod,qmbasis=qmbasis,addkw=addkw,qmmem=qmmem)
    with open(logfile,'a') as f:
        print("Total S",sum(Smol.values())-len(Smol)+1,file=f)
    g16tb=printtb(nproc,d,atomtype,atomxyz,Smol,r=r,qmmethod=qmmethod,qmbasis=qmbasis,addkw=addkw,qmmem=qmmem)
    return g16mol+g16tb

def readforce(jobname):
    d={}
    with open(jobname+".id") as f:
        for line in f:
            s=line.split()
            d[s[0]]=s[1]
    with open(jobname+".log") as f:
        b =False
        atoms={}
        for line in f:
            if "Forces (Hartrees/Bohr)" in line:
                b=True
                i=0
                continue
            if b:
                i+=1
                if i>2:
                    if line.startswith(" ------------"):
                        break
                    s=line.split()
                    #print(d[s[0]],s[2],s[3],s[4])
                    atoms[int(d[s[0]])]=np.array([float(x) for x in s[2:5]])
    return atoms

def takeforce(g16,logfile):
    atomforce={}
    i=0
    while "g16 mol"+str(i+1)+".gjf" in g16:
        i+=1
        forces=readforce("mol"+str(i))
        for atom,force in forces.items():
            atomforce[atom]=force
    n=i
    twobodyforce={}
    for i in range(1,n+1):
        for j in range(i+1,n+1):
            if not ("g16 tb"+str(i)+"-"+str(j)+".gjf" in g16):
                continue
            forces=readforce("tb"+str(i)+"-"+str(j))
            for atom,force in forces.items():
                twobodyforce[atom]=twobodyforce[atom]+force-atomforce[atom] if atom in twobodyforce else force-atomforce[atom]
    finalforces=[]
    for atom,force in sorted(atomforce.items(),key=lambda item:item[0]):
        finalforce=force+twobodyforce[atom] if atom in twobodyforce else force
        print ("".join("%16.9f"%x for x in finalforce))
        finalforces.append(finalforce)
    with open(logfile,'a') as f:
        forcesum=np.sum(finalforces,axis=0)
        print("".join("%16.9f"%x for x in forcesum),file=f,end='')
        forcesumdis=np.sqrt(np.sum(np.square(forcesum)))
        print("%16.9f"%forcesumdis,file=f)
        print("\n",file=f)

def run(nproc_sum,nproc,r,xyzfilename,pdbfilename,qmmethod,qmbasis,addkw,qmmem,atombondnumber,logfile):                
    os.system('obabel -ixyz '+xyzfilename+' -opdb -O '+pdbfilename+' > /dev/null')
    g16commands=readbond(nproc=nproc,pdbfilename=pdbfilename,r=r,qmmethod=qmmethod,qmbasis=qmbasis,addkw=addkw,qmmem=qmmem,atombondnumber=atombondnumber,logfile=logfile)
    n=nproc_sum//nproc    
    os.system('Nproc='+str(n)+';Pfifo=/tmp/$$.fifo;mkfifo $Pfifo;exec 6<>$Pfifo;rm -f $Pfifo;for ((i=1;i<=$Nproc;i++));do echo;done >&6;'+''.join(['read -u6 ; { '+g16command +' ; echo >&6; }&' for g16command in g16commands])+'wait;exec 6>&-')
    takeforce(g16=g16commands,logfile=logfile) 

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
    run(nproc_sum=int(sys.argv[1]),nproc=qmproc,r=dl,xyzfilename=xyzfilename,pdbfilename=pdbfilename,qmmethod=qmmethod,qmbasis=qmbasis,addkw=addkw,qmmem=qmmem,atombondnumber=atombondnumber,logfile="force.log")
