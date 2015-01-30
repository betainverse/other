#!/usr/bin/python
"""
countSeats.py
Given: length of a row of subway seats

Return: Mean and standard deviation of the fraction of occupied seats, assuming that people choose seats randomly but are unwilling to sit next to another person, until no acceptable seats remain.
"""
from random import randint
from numpy import mean, std
from datetime import datetime

total=0
available = []
occupied = []

# Given a list of occupied seats, and number of total seats,
# return a list of available seats.
# getAvailable([0,2],6) -> [4,5]
#def getAvailable(occupied,total):
#    available = []
#    for i in range(total):
#        if i not in occupied and i-1 not in occupied and i+1 not in occupied:
#            available.append(i)
#    return available

# Given a list of occupied seats, and a number of total seats,
# Add a person, and return a new list of occupied seats.
# Assume available seats remain
# addPerson([0,2],6) -> [0,2,4] or [0,2,5]
#def addPerson(occupied,total):
#    available = getAvailable(occupied,total)
#    remaining = len(available)
#    occupied.append(available[randint(0,remaining-1)])
#    return occupied
def addPerson():
    global total,occupied,available
    remaining = len(available)
    chosenSeat = available[randint(0,remaining-1)]
    occupied.append(chosenSeat)
    available.remove(chosenSeat)
    if chosenSeat+1 in available:
        available.remove(chosenSeat+1)
    if chosenSeat-1 in available:
        available.remove(chosenSeat-1)

# Given a total number of seats, return a list of seats that become occupied.
# fillSeats(6) -> one of [0,2,4],[0,2,5],[1,3,5],etc.
#def fillSeats(total):
#    occupied = []
#    available = getAvailable(occupied,total)
#    while len(available)>0:
#        occupied = addPerson(occupied,total)
#        available = getAvailable(occupied,total)
#    return occupied

def fillSeats():
    global total, available, occupied
    while len(available)>0:
        addPerson()

def emptySeats():
    global total, available, occupied
    occupied = []
    available = range(total)

# Given a total number of seats S and a number of tries n,
# try filling S seats at random n times, and return a list of the number of
# occupied seats for all tries.
# def collectTries(S,n):
#     trials = []
#     for i in range(n):
#         occupied = fillSeats(S)
#         trials.append(len(occupied))
#     return trials

# def collectTriesDict(S,n):
#     trials = {}
#     for i in range(n):
#         occupied = len(fillSeats(S))
#         if occupied in trials.keys():
#             trials[occupied] = trials[occupied]+1
#         else:
#             trials[occupied] = 1
#     return trials

def collectTries(n):
    global total, available, occupied
    trials = []
    for i in range(n):
        emptySeats()
        fillSeats()
        trials.append(float(len(occupied))/total)
    return mean(trials),std(trials)

# Given a total number of seats S and a number of tries n,
# calculate the mean and standard deviation for the fraction of
# occupied seats.

# def getMeanStd(S,n):
#     fracOccupied = [float(x)/S for x in collectTries(S,n)]
#     return mean(fracOccupied),std(fracOccupied)

# def getMeanStdDict(totalSeats,totalObservations):
#     count = collectTriesDict(totalSeats,totalObservations)
#     meanFrac = sum([float(seats)*count[seats]/totalSeats for seats in count.keys()])/totalObservations
#     stdFrac = sum([count[seats]*abs(meanFrac-float(seats)/totalSeats) for seats in count.keys()])/totalObservations
#     return meanFrac,stdFrac

def main():
    
    #print getAvailable([0,2],6)
    #print addPerson([0,2],6)
    #print fillSeats(6)
    #print collectTries(25,50)
    #print getMeanStd(25,100)
    #print getMeanStd(25,1000)
    #print getMeanStd(25,10000)
    #print getMeanStd(25,100000)
    #print getMeanStd(25,1000000)
    #print getMeanStdDict(25,1000)
    #print getMeanStdDict(25,10000)
    #print getMeanStdDict(25,100000)
    #print getMeanStdDict(25,1000000)
    global total, available, occupied
    total=25
    # print datetime.now().time()
    # print 1000, collectTries(10000)
    # print datetime.now().time()
    # print 10000, collectTries(100000)
    # print datetime.now().time()
    # print 100000, collectTries(1000000)
    # print datetime.now().time()
    # print 1000000, collectTries(10000000)
    print datetime.now().time()
    print 10000000, collectTries(100000000)
    print datetime.now().time()
main()
