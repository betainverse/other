#!/usr/bin/python
"""
PREdist.py computes the distances between the CB of PRE-labeled residues and the amide N of the binding partner.
"""

from sys import stdin
import math

global mutant_residues 
mutant_residues = [78,114,139,145,167,189,194,205,218]
mutant_residues.sort()

class Atom:
    def __init__(self, PDBline):
        fields = PDBline.split()
        self.atnum = int(fields[1])
        self.atname = fields[2]
        self.resname = fields[3]
        self.resnum = int(fields[4])
        self.x = float(fields[5])
        self.y = float(fields[6])
        self.z = float(fields[7])
    def dist(self,other):
        return math.sqrt((self.x-other.x)**2+(self.y-other.y)**2+(self.z-other.z)**2)

def parse_atoms(openfile):
    amides = {}
    mutants = {}
    for line in openfile.readlines():
        if 'ATOM' in line:
            newatom = Atom(line)
            if newatom.atname == 'N' and newatom.resnum >= 1238:
                amides[newatom.resnum]=newatom
            elif newatom.atname == 'CB' and newatom.resnum in mutant_residues:
                mutants[newatom.resnum]=newatom
    return amides,mutants

def print_distances(amides,mutants):
    amide_residues = amides.keys()
    amide_residues.sort()
    header = ""
    for residue in mutant_residues:
        header+= "\t%d"%residue
    print header
    for res1 in amide_residues:
        distance_string = str(res1)
        for res2 in mutant_residues:
            distance_string+="\t%0.3f"%mutants[res2].dist(amides[res1])
        print distance_string
        

def main():
    amides,mutants=parse_atoms(stdin)
    print_distances(amides,mutants)
                
            

if __name__ == "__main__":
    main()
