#!/usr/bin/env python3
# Add root of readbond.py to PATH
# Usage:
# readbond.py 28 xxx.xyz
from AIMDBlock import AIMDFragmentation
import time
import logging
import sys

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
    time1=time.time()
    xyzfilename=sys.argv[2]
    pdbfilename=xyzfilename+".pdb"
    atombondnumber={"C":4,"H":1,"O":2}
    qmproc,qmmethod,qmbasis,addkw,dl,qmmem=readmfccin(mfccinfilename="mfcc.in")
    AIMDFragmentation(nproc_sum=int(sys.argv[1]),nproc=qmproc,cutoff=dl,xyzfilename=xyzfilename,pdbfilename=pdbfilename,qmmethod=qmmethod,qmbasis=qmbasis,addkw=addkw,qmmem=qmmem,atombondnumber=atombondnumber,logfile="force.log").run()
    time2=time.time()
    logging.warning("Time consumed: "+str(time2-time1))
