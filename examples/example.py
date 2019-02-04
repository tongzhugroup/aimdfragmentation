from aimdfragmentation import AIMDFragmentation
AIMDFragmentation(nproc_sum=28,
                  nproc=4,
                  cutoff=4.0,
                  xyzfilename="ch4.xyz",
                  qmmethod="mn15",
                  qmbasis="6-31g(d)",
                  qmmem="400MW",
                  ).run()
