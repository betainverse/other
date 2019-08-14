#!/usr/bin/env python
import csv
import matplotlib.pyplot as plt
import numpy as np

# user-editable fields. See syringe program below if necessary.
upper_anisotropy_limit = 0.163
lower_anisotropy_limit = 0.141
baseline_start_time = 0
injection_start_time = 1201
data_end_time = 100000
cuvette_volume_uL = 3015
protein_syringe_conc_uM = 3
syringe_program_uL = [1,1,1,1,1,1,2,2,2,2,2,2,5,5,5,5,5,5,10,10,10,10,10]
#syringe_program_uL = [10,10,10,10,10,10,10,10,10,10]
time_between_injections_s = 240
exclude_times = [] # this has not been implemented - need to address sig figs

def main():
    openfile = open('20190807 reduced C9R tris 25C 750mM run 1.ifx')
    lines = openfile.readlines()
    start_data = get_start_data(lines)
    title,comment = get_title_comment(lines)
    title=title[0:-2]
    time2anisotropy = line2timeanisotropy(lines[start_data:])
    plot_time_anisotropy(time2anisotropy,title)
    #print start_data
    #time2volume = time2injected_volume(time2anisotropy)
    times = time2anisotropy.keys()
    conc2aniso = concentration2anisotropies(time2anisotropy)
    plot_conc_anisotropy(conc2aniso,title)
    plot_conc_anisotropy_avg_std(conc2aniso,title)
    write_dynafit_input(conc2aniso,title,comment)

def get_title_comment(lines):
    for line in lines:
        if "Title" in line:
            title = line.split('=')[1]
        if "Comment" in line:
            comment = line.split('=')[1]
            break
    return (title,comment)
def get_start_data(lines):
    for i in range(len(lines)):
        if "Columns" in lines[i]:
            columns= lines[i].split('=')[1]
            fields = columns.split(',')
            if 'Time' not in fields[0]:
                print "error: time field not located"
            elif 'Anisotropy' not in fields[2]:
                print "error: anisotropy field not located"                
        elif "[Data]" in lines[i]:
            start_data=i+1
            break
    return start_data

# given a list of lines of data only, with headers removed, return a dictionary of times and anisotropies
# automatically remove and report egregiously bad data.
def line2timeanisotropy(lines):
    time2anisotropy = {}
    for line in lines:
        fields = line.split('\t')
        time = float(fields[0])
        anisotropy = float(fields[2])
        if time > baseline_start_time and time < data_end_time and time not in exclude_times:
            if anisotropy > lower_anisotropy_limit and anisotropy < upper_anisotropy_limit:
                time2anisotropy[time]=anisotropy
            else:
                print "Warning: Anisotropy %0.4f at time %0.2f is outside standard bounds and has been excluded."%(anisotropy,time)
    return time2anisotropy

def plot_time_anisotropy(time2anisotropy,title):
    times = time2anisotropy.keys()
    times.sort()
    anisotropies = [time2anisotropy[time] for time in times]
    fig = plt.figure()
    plt.plot(times,anisotropies)
    plt.xlabel('Time (s)')
    plt.ylabel('Anisotropy')
    plt.title(title)
    plt.show()
    fig.savefig(title+' time_v_anisotropy.pdf',bbox_inches='tight')

def time2injected_volume(times2anisotropy):
    #syringe_program_uL = [1,1,1,1,1,1,2,2,2,2,2,2,5,5,5,5,5,5,10,10,10,10,10]
    times = times2anisotropy.keys()
    times.sort()
    time2volume = {}
    injection_times = [injection_start_time+time_between_injections_s*x for x in range(len(syringe_program_uL))]
    for time in times:
        if time < injection_times[0]:
            time2volume[time]=0
        elif time > injection_times[-1]:
            time2volume[time]=sum(syringe_program_uL)
        else:
            for i in range(1,len(injection_times)):
                if time > injection_times[i-1] and time < injection_times[i]:
                    time2volume[time]=sum(syringe_program_uL[0:i])
                    break
    return time2volume

def time2concentration(time2volume):
    time2conc={}
    times = time2volume.keys()
    times.sort()
    for time in times:
        injected_volume = time2volume[time]
        total_volume = float(injected_volume + cuvette_volume_uL)
        numerator = float(injected_volume * protein_syringe_conc_uM)
        concentration_nM = 1000*numerator/total_volume
        time2conc[time]=concentration_nM
    return time2conc

def concentration2anisotropies(time2anisotropy):
    times = time2anisotropy.keys()
    times.sort()
    time2conc = time2concentration(time2injected_volume(time2anisotropy))
    #print time2conc
    newtimes = time2conc.keys()
    newtimes.sort()
    #print times
    concentration2anisotropies = {}
    for time in times:
        conc = time2conc[time]
        if conc in concentration2anisotropies.keys():
            concentration2anisotropies[conc].append(time2anisotropy[time])
        else:
            concentration2anisotropies[conc]=[time2anisotropy[time]]
    return concentration2anisotropies

def plot_conc_anisotropy(concentration2anisotropies,title):
    concentrations = concentration2anisotropies.keys()
    concentrations.sort()
    xlist = []
    ylist = []
    for conc in concentrations:
        for anisotropy in concentration2anisotropies[conc]:
            xlist.append(conc)
            ylist.append(anisotropy)
    fig = plt.figure()
    plt.scatter(xlist,ylist)
    plt.xlabel('Tetramer concentration (nM)')
    plt.ylabel('Anisotropy')
    plt.title(title)
    plt.show()
    fig.savefig(title+' conc_v_anisotropies.pdf',bbox_inches='tight')

def plot_conc_anisotropy_avg_std(concentration2anisotropies,title):
    concentrations = concentration2anisotropies.keys()
    concentrations.sort()
    means = []
    stds = []
    for conc in concentrations:
        anisotropy = np.mean(concentration2anisotropies[conc])
        stdev = np.std(concentration2anisotropies[conc])
        means.append(anisotropy)
        stds.append(stdev)
    fig = plt.figure()
    plt.errorbar(concentrations,means,yerr=stds,fmt='o')
    plt.xlabel('Tetramer concentration (nM)')
    plt.ylabel('Anisotropy')
    plt.title(title)
    plt.show()
    fig.savefig(title+' conc_v_anisotropy_mean_std.pdf',bbox_inches='tight')
        

def write_dynafit_input(concentration2anisotropies,title,comment):
    concentrations = concentration2anisotropies.keys()
    concentrations.sort()
    means = []
    stds = []
    openfile = open(title+'_dynafit_input.txt','w')
    openfile.write(';'+comment) #comment has the windows newline char
    openfile.write('Concentration\tAverage Anisotropy\tStDev Anisotropy\r\n')
    for conc in concentrations:
        anisotropy = np.mean(concentration2anisotropies[conc])
        stdev = np.std(concentration2anisotropies[conc])
        means.append(anisotropy)
        stds.append(stdev)
        openfile.write('%0.8f\t%0.8f\t%0.8f\r\n'%(conc,anisotropy,stdev))
    minimum = min(means)
    maximum = max(means)
    midpoint = (minimum+maximum)/2
    deltar = maximum-minimum
    openfile.write('Minimum: %0.6f\r\nMidpoint: %0.6f\r\nMaximum: %0.6f\r\nDelta_r: %0.6f\r\n'%(minimum,midpoint,maximum,deltar))
    openfile.close()
                    
    

main()

