#!/usr/bin/python
"""
PREdist.py computes the distances between the CB of PRE-labeled residues and the amide N of the binding partner.
"""

from sys import stdin,argv
import math
import pylab
import matplotlib.lines as mpllines
import matplotlib.ticker as mplticker


global mutant_residues 
mutant_residues = [78,114,139,145,167,189,194,205,218]
#mutant_residues = [78,139,145,218]
#mutant_residues = [114,167,189,194]
mutant_residues.sort()
global fignum
fignum = 1

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

#def plot_distances(amides,mutant_CB,outfile):
#    amide_residues = amides.keys()
#    amide_residues.sort()
#    distances = [amides[res].dist(mutant_CB) for res in amide_residues]
#    thresh_distances = [dist if dist<30 else 30 for dist in distances]
#    mybargraph(amide_residues,thresh_distances,outfile)

def get_distances(amides,mutant_CB):
    amide_residues = amides.keys()
    amide_residues.sort()
    distances = [amides[res].dist(mutant_CB) for res in amide_residues]
    return distances

def mybargraph(xlist,ylist,outfile):
    fig = pylab.figure(figsize=(10,3))
    pylab.bar(xlist,ylist)
    #pylab.plot(xlist,thresh,'r')
    #pylab.plot(xlist,thresh2,'r')
    pylab.axis(xmin=min(xlist),xmax=max(xlist))
    pylab.xlabel("Residue")
    pylab.ylabel("Distance")
    ax = fig.gca()
    ax.xaxis.set_major_locator(mplticker.MultipleLocator(10))
    for line in ax.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    pylab.gcf().subplots_adjust(bottom=0.15)
    pylab.savefig(outfile)

def plot_all_distances(amides,mutants,outfile):
    amide_residues = amides.keys()
    amide_residues.sort()
    fig = pylab.figure(1)
    numplots = len(mutant_residues)
    for i in range(numplots):
        print mutant_residues[i]
        subplotidentifier = str(numplots)+str(1)+str(i+1)
        pylab.subplot(int(subplotidentifier))
        singlebargraph(amide_residues,get_distances(amides,mutants[mutant_residues[i]]),mutant_residues[i],fig)
    pylab.savefig(outfile,orientation='portrait')

def plot_some_distances(amides,mutants,mutant_subset,outfile):
    global fignum
    amide_residues = amides.keys()
    amide_residues.sort()
    fig = pylab.figure(fignum)
    fig.suptitle(outfile)
    fignum+=1
    numplots = len(mutant_subset)
    for i in range(numplots):
        print mutant_subset[i]
        subplotidentifier = str(numplots)+str(1)+str(i+1)
        pylab.subplot(int(subplotidentifier))
        singlebargraph(amide_residues,get_distances(amides,mutants[mutant_subset[i]]),mutant_subset[i],fig)
    pylab.savefig(outfile)


def singlebargraph(xlist,ylist,name,figure):
    fig = figure
    pylab.bar(xlist,ylist,color='w')
    #pylab.plot(xlist,thresh,'r')
    #pylab.plot(xlist,thresh2,'r')
    pylab.axis(xmin=min(xlist),xmax=max(xlist))
    pylab.axis(ymin=0,ymax=30)
    pylab.xlabel("Residue")
    pylab.ylabel("Dist "+str(name))
    ax = fig.gca()
    ax.xaxis.set_major_locator(mplticker.MultipleLocator(10))
    ax.yaxis.set_major_locator(mplticker.MultipleLocator(10))
    for line in ax.get_xticklines():
        line.set_marker(mpllines.TICKDOWN)
    labels = ax.get_xticklabels()
    for label in labels:
        label.set_rotation(30) 
    #pylab.gcf().subplots_adjust(bottom=0.15)
    #pylab.savefig(outfile)
    

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
    mutant_subset1 = [78,139,145,218]
    mutant_subset2 = [114,167,189,194]
    #infile = argv[1]
    for infile in argv[1:]:
        openfile = open(infile,'r')
        amides,mutants=parse_atoms(openfile)
        openfile.close()
        plot_some_distances(amides,mutants,mutant_subset1,infile.split('.')[0]+'_1.pdf')
        plot_some_distances(amides,mutants,mutant_subset2,infile.split('.')[0]+'_2.pdf')
            

if __name__ == "__main__":
    main()
