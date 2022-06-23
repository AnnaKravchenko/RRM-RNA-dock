#!/usr/bin/python2.7
# This script is based on the Michael J. Harms's script of the same name 

import os, sys

def pdbAtomRenumber(pdb):
    """
    Renumber all atoms in pdb file, starting from 1.
    """

    out = []
    counter = 1
    for line in pdb:
        # For and ATOM record, update residue number
        if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
            out.append("%s%5s%s" % (line[0:6],counter,line[11:]))
            counter += 1
        else:
            out.append(line)
    return out

pdb_file = sys.argv[1]
new_name = sys.argv[2]
   
# Read in the pdb file
f = open(pdb_file,'r')
pdb = f.readlines()
f.close()

out = pdbAtomRenumber(pdb)

out_file = new_name
g = open(out_file,'w')
g.writelines(out)
g.close()