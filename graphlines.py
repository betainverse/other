#!/usr/bin/env python
"""
Read in a csv file with wavelengths in the first column, and absorption
values in the remaining columns, with a header line giving molar equivalents
of additive for each column. 
Output overlaid spectra.
Possibly do some baseline correction.
Possibly extract a particular wavelength
"""
import sys
from pylab import *

def subtractbaseline(listofabs):
    subspectrum = listofabs[510:620]
    newzero = min(subspectrum)-0.01
    return [x-newzero for x in listofabs]

def main():
    if len(sys.argv) != 3:
        print "Usage:"
        print "graphlines.py spectra.csv chart.pdf"
        return
    datadict = {}
    infile = sys.argv[1]
    outfile = sys.argv[2]
    openfile = open(infile,'r')
    header = openfile.readline()
    equivalents = header[:-1].split(',')[1:]
    print equivalents
    for line in openfile.readlines():
        columns = line[:-1].split(',')
        datadict[int(columns[0])]=columns[1:]
    openfile.close()
    spectradict={}
    wavelengths = datadict.keys()
    wavelengths.sort()
    for i in range(len(equivalents)):
        key = float(equivalents[i])
        spectradict[key]=[]
        for wavelength in wavelengths:
            absorbance = float(datadict[wavelength][i])
            spectradict[key].append(absorbance)
    correctedbaselinedict = {}
    for key in spectradict.keys():
        correctedbaselinedict[key]=subtractbaseline(spectradict[key])
    fig = figure()
    for eq in spectradict.keys():
        plot(wavelengths,correctedbaselinedict[eq])
    ax = fig.gca()
    ax.set_xlim(300,800)
    ax.set_ylim(0,0.12)
    savefig(outfile)
            
    

main()
