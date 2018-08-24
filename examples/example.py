from AIMDBlock import AIMDBlock
AIMDBlock(  nproc_sum=28,
            nproc=4,
            cutoff=2.0,
            xyzfilename="ch4.xyz",
            pdbfilename="ch4.pdb",
            qmmethod="mn15",
            qmbasis="6-31g(d)",
            addkw="scf=xqc",
            qmmem="400MW",
            atombondnumber={"C":4,"H":1,"O":2},
            logfile="force.log"
        ).run()
