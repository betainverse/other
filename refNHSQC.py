#!/usr/bin/python
"""
ref.py takes a reference value, a procpar file, and a fid.com file,
and rewrites the fid.com file with referenced values.

Usage: ref.py 599.7981547 procpar fid.com
"""

import sys
import os

if len(sys.argv) not in [2,4]:
    print "Usage: ref.py 599.7981547 procpar fid.com\n"
    exit
elif len(sys.argv) == 2:
    print "\nUsing existing procpar and fid.com files by default.\n"
    ppfile = open('procpar','r')
    fidfile = open('fid.com','r')
elif len(sys.argv) == 4:
    ppfile = open(sys.argv[2],'r')
    fidfile = open(sys.argv[3],'r')

h0 = float(sys.argv[1])
c0 = 0.251449530 * h0
n0 = 0.101329118 * h0
p0 = 0.404808636 * h0

fidlines = fidfile.readlines()
fidfile.close()
pplines = ppfile.readlines()
ppfile.close()

def getparam(param):
    for i in range(len(pplines)):
        if pplines[i].split()[0] == param:
            return float(pplines[i+1].split()[1])



H_carrier_Hz = getparam('sfrq')
H_carrier_ppm = 1000000 * ( H_carrier_Hz - h0 )/h0
tof = getparam('tof')
print "PROTON:\ntof  = %8.2f \t carrier = %11.7f ppm = %11.7f MHz"%(tof,H_carrier_ppm,H_carrier_Hz)
C_carrier_Hz = getparam('dfrq')
C_carrier_ppm = 1000000 * ( C_carrier_Hz - c0 )/c0
dof = getparam('dof')
print "CARBON:\ndof  = %8.2f \t carrier = %11.7f ppm = %11.7f MHz"%(dof,C_carrier_ppm,C_carrier_Hz)
N_carrier_Hz = getparam('dfrq2')
N_carrier_ppm = 1000000 * ( N_carrier_Hz - n0 )/n0
dof2 = getparam('dof2')
print "NITROGEN:\ndof2 = %8.2f \t carrier = %11.7f ppm = %11.7f MHz"%(dof2,N_carrier_ppm,N_carrier_Hz)

fidfile = open('newfid.com','w')

for line in fidlines:
    columns = line.split()
    if len(columns) < 4:
        fidfile.write(line)
    elif columns[0] == '-xOBS':
        columns[1] = '%11.7f'%H_carrier_Hz
        columns[3] = '%11.7f'%N_carrier_Hz
        newline = '  '+'    '.join(columns)+'\n'
        fidfile.write(newline)
    elif columns[0] == '-xCAR':
        columns[1] = '%11.7f'%H_carrier_ppm
        columns[3] = '%11.7f'%N_carrier_ppm
        newline = '  '+'    '.join(columns)+'\n'
        fidfile.write(newline)
    else:
        fidfile.write(line)
fidfile.write('#Referenced using spcfrq = %11.7f and refNHSQC.py\n\n'%h0)
fidfile.close()    

#st = os.stat('newfid.com')
#os.chmod('newfid.com', st.st_mode | 0111)
os.system('chmod u+x newfid.com')
os.system('./newfid.com')
