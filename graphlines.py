#!/usr/bin/env python
"""
Read in a set of csv files from the UV-vis spectrophotometer.
Output a pdf with overlaid spectra.
Possibly do some baseline correction.
"""
import sys
from pylab import *

outputfile = 'SczACoTitration2.pdf'
inputdirectory = '/home/edmonds/Documents/20131209SczACo/20131209CoSczA/'
# inputfiles = ['SCZA0R1.CSV',
#               'SCZA1R1.CSV',
#               'SCZA2R1.CSV',
#               'SCZA3R1.CSV',
#               'SCZA4R1.CSV',
#               'SCZA5R1.CSV',
#               'SCZA6R1.CSV',
#               'SCZA7R1.CSV',
#               'SCZA8R1.CSV',
#               'SCZA9R1.CSV',
#               'SCZA10R1.CSV',
#               'SCZA11R1.CSV',
#               'SCZA12R1.CSV',
#               'SCZA13R1.CSV',
#               'SCZA14R1.CSV',
#               'SCZA15R1.CSV',
#               'SCZA16R1.CSV',
#               'SCZA17R1.CSV',
#               'SCZA18R1.CSV']
# spectrafiles is a dictionary relating number of molar equivalents of 
# added cobalt to the filename
spectrafiles = {0.0:'SCZA0R1.CSV',
                0.145:'SCZA1R1.CSV',
                0.29:'SCZA2R1.CSV',
                0.435:'SCZA3R1.CSV',
                0.58:'SCZA4R1.CSV',
                0.725:'SCZA5R1.CSV',
                0.87:'SCZA6R1.CSV',
                1.015:'SCZA7R1.CSV',
                1.16:'SCZA8R1.CSV',
                1.305:'SCZA9R1.CSV',
                1.45:'SCZA10R1.CSV',
                1.595:'SCZA11R1.CSV',
                1.74:'SCZA12R1.CSV',
                1.885:'SCZA13R1.CSV',
                2.03:'SCZA14R1.CSV',
                2.175:'SCZA15R1.CSV',
                2.32:'SCZA16R1.CSV',
                2.68:'SCZA17R1.CSV',
                3.045:'SCZA18R1.CSV'}

proteinconc = 0.000255 #255uM = 0.255mM = 0.000255M

def subtractbaseline(listofabs):
    subspectrum = listofabs[510:620]
    newzero = min(subspectrum)-0.01
    return [x-newzero for x in listofabs]

def getXaxis(filename):
    Xvalues = []
    openfile = open(filename,'r')
    header = openfile.readline()
    for line in openfile.readlines():
        columns = line[:-1].split(',')
        x = int(columns[0])
        Xvalues.append(x)
    openfile.close()
    return Xvalues

def getYaxis(filename,proteinconc):
    Yvalues = []
    openfile = open(filename,'r')
    header = openfile.readline()
    for line in openfile.readlines():
        columns = line[:-1].split(',')
        y = float(columns[1])/proteinconc
        Yvalues.append(y)
    openfile.close()
    return Yvalues

def getYval(wavelength,yValues,xValues):
    idx = xValues.index(wavelength)
    return yValues[idx]

def plotTitration(curvedict,xValues,wavelength,outputfile):
    equivalents = curvedict.keys()
    equivalents.sort()
    xlist = []
    ylist = []
    for eq in equivalents:
        xlist.append(eq)
        ylist.append(getYval(wavelength,curvedict[eq],xValues))
    fig = figure()
    plot(xlist,ylist)
    plot(xlist,ylist,'bo')
    ylabel("$\epsilon_{%d}$ $(M^{-1}cm^{-1})$"%wavelength)
    xlabel("Molar equivalents of Co added to SczA")
    ax = fig.gca()
    ax.set_xlim(0,3.1)
    ax.set_ylim(-10,350)
    basefilename = outputfile.split('.')[0]
    savefig(basefilename+'-%d.pdf'%wavelength)
     

def main():
    if len(sys.argv) != 1:
        print "Usage:"
        print "Enter your information in script header, then run:"
        print "graphlines.py"
        return
    inputfiles = spectrafiles.values()
    xValues = getXaxis(inputdirectory+inputfiles[0])
    #curves = []
    curvedict = {}
    equivalents = spectrafiles.keys()
    equivalents.sort()
    for eq in equivalents:
        filename = spectrafiles[eq]
        yValues = getYaxis(inputdirectory+filename,proteinconc)
        #curves.append(yValues)
        curvedict[eq]=yValues
    fig = figure()
#    for curve in curves:
#        plot(xValues,curve)
    for eq in equivalents:
        plot(xValues,curvedict[eq])
    ax = fig.gca()
    ax.set_xlim(300,800)
    ax.set_ylim(-10,400)
    ylabel("$\epsilon$ $(M^{-1}cm^{-1})$")
    xlabel("$\lambda$ (nm)")
    savefig(outputfile)
    plotTitration(curvedict,xValues,560,outputfile)        
    

main()
