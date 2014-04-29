#!/usr/bin/python
"""
modifyPrintSparky.py takes a .ps printout of a spectrum/overlay from sparky and modifies several parameters to make nicer output.
"""

from sys import stdin,argv

def setFalse(line):
    pos = line.find('true')
    return line[:pos]+'false'+line[pos+4:]

def doublethickness(line):
    pos = line.find('1.0')
    return line[:pos]+'2.0'+line[pos+2:]

def main():
    infile = argv[1]
    openfile = open(infile,'r')
    for line in openfile.readlines():
        newline = line
        if 'true' in line:
            if 'topAxisScale' in line or 'rightAxisScale' in line or 'topTick' in line or 'rightTick' in line:
                newline = setFalse(line)
        elif '1.0 setlinewidth' in line:
            newline = doublethickness(line)
        print newline[:-1]
    openfile.close()

if __name__ == "__main__":
    main()
