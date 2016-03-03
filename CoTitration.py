#!/usr/bin/env python
"""
Read in a set of csv files from the UV-vis spectrophotometer.
Output a pdf with overlaid spectra.
Possibly do some baseline correction.
"""
import sys
from pylab import *
from scipy.optimize import curve_fit
import numpy as np

outputname = 'SczAK9ACoTitration'
inputdirectory = '/Users/edmonds/Dropbox/giedroc lab/SczA/Co titrations/20160229/K9A/'

# spectrafiles is a dictionary relating number of injections of 
# added cobalt to the filename
spectrafiles = {0.0:'1.CSV',
                1:'2.CSV',
                2:'3.CSV',
                3:'4.CSV',
                4:'5.CSV',
                6:'6.CSV'}

initial_protein_conc = 0.000372 #318uM
initial_volume = 150 #uL
injection_volume = 1.79  #uL
Co_stock_conc = 0.01506 # 15060uM
#vol_per_Co_equiv = 1.79*2 #uL
#vol_per_Co_equiv = initial_protein_conc*initial_volume/Co_stock_conc
#print vol_per_Co_equiv, "uL"

outputpdf = outputname + '.pdf'
outputsheet = outputname + '.csv'
openfile = open(outputsheet,'w')

def getXvalues(filename):
    Xvalues = []
    openfile = open(filename,'r')
    header = openfile.readline()
    for line in openfile.readlines():
        columns = line[:-1].split(',')
        x = int(columns[0])
        Xvalues.append(x)
    openfile.close()
    return Xvalues

def getYvalues(filename):
    Yvalues = []
    openfile = open(filename,'r')
    header = openfile.readline()
    for line in openfile.readlines():
        columns = line[:-1].split(',')
        y = float(columns[1])
        Yvalues.append(y)
    openfile.close()
    return Yvalues

def combineXY(Xvalues,Yvalues):
    return dict(zip(Xvalues,Yvalues))

def subtractBaseline(XYdict):
    baseline = XYdict[800]
    newdict = {}
    for x in XYdict.keys():
        newdict[x] = XYdict[x]-baseline
    return newdict

def fitfunc(x,a,b,c):
    return a*np.exp(-b*x)+c

#def fitfunc(x,a,b,c):
#    return -a*np.log10(1/(1-(b*(x**-4))))+c

def correctBaseline(XYdict,filename):
    xValues = [x for x in XYdict.keys() if x in range(315,425) or x in range(680,900)]
    xValues.sort()
    yValues = [XYdict[x] for x in xValues]
    popt, pcov = curve_fit(fitfunc,xValues,yValues,[0,0,0],maxfev=10000)
    print popt
    xx = np.linspace(300,900,50)
    yy = [fitfunc(x,*popt) for x in xx]
    fig=figure()
    #plot(xValues,yValues)
    plot(XYdict.keys(),XYdict.values())
    plot(xx,yy,'r')
    ax = fig.gca()
    ax.set_xlim(300,900)
    ax.set_ylim(-0.1,0.1)
    savefig(filename+'.pdf')
    newdict = {}
    for x in XYdict.keys():
        newdict[x] = XYdict[x]-fitfunc(x,*popt)
    return newdict
    

def convertByMolarity(XYdict,injections):
    final_vol = initial_volume+injections*injection_volume
    protein_conc = initial_protein_conc*initial_volume/final_vol
    newdict = {}
    for x in XYdict.keys():
        newdict[x] = XYdict[x]/protein_conc
    return newdict

def importfile(injections,filename):
    rawData = combineXY(getXvalues(inputdirectory+filename),getYvalues(inputdirectory+filename))
    #return convertByMolarity(subtractBaseline(rawData),injections)
    return convertByMolarity(correctBaseline(rawData,filename),injections)
    #return convertByMolarity(rawData,injections)

def printHeader(spectrafiles):
    injections = spectrafiles.keys()
    injections.sort()
    final_volumes = [initial_volume + x*injection_volume for x in injections]
    protein_concentrations = [initial_protein_conc*initial_volume/v for v in final_volumes]
    co_concentrations = [x*injection_volume*Co_stock_conc/(initial_volume + x*injection_volume) for x in injections]
    firstHeader = 'injections Co added,'+','.join(['%d'%inj for inj in injections])+'\n'
    secondHeader = 'protein conc,'+','.join(['%0.1f'%(1000000*conc) for conc in protein_concentrations])+'\n'
    thirdHeader = 'Co conc,'+','.join(['%0.1f'%(1000000*conc) for conc in co_concentrations])+'\n'
    fourthHeader = 'wavelength/volume,' +','.join(['%0.3f'%v for v in final_volumes])+'\n'
    #print firstHeader
    #print secondHeader
    openfile.write(firstHeader)
    openfile.write(secondHeader)
    openfile.write(thirdHeader)
    openfile.write(fourthHeader)

def printline(wavelength,allData):
    injections = allData.keys()
    injections.sort()
    abs_values = [allData[key][wavelength] for key in injections]
    nextLine = '%d,'%wavelength+','.join(['%0.8f'%a for a in abs_values])+'\n'
    #print nextLine
    openfile.write(nextLine)
    
def plotSpectra(allData):
    fig = figure()
    xValues = allData[0].keys()
    xValues.sort()
    keys = allData.keys()
    keys.sort()
    for key in keys:
        plot(xValues,[allData[key][x] for x in xValues])
    ax = fig.gca()
    ax.set_xlim(320,800)
    ax.set_ylim(-10,350)
    ylabel("$\epsilon$ $(M^{-1}cm^{-1})$")
    xlabel("$\lambda$ (nm)")
    savefig(outputpdf)
    

def main():
    if len(sys.argv) != 1:
        print "Usage:"
        print "Enter your information in script header, then run:"
        print "graphlines.py"
        return
    allData = {}
    for key in spectrafiles.keys():
        filename = spectrafiles[key]
        allData[key]=importfile(key,filename)
    printHeader(spectrafiles)
    wavelengths = allData[0].keys()
    wavelengths.sort()
    for wavelength in wavelengths:
        printline(wavelength,allData)
    plotSpectra(allData)
    openfile.close()

main()
