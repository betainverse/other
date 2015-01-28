#!/usr/bin/python
"""
countSeats.py
Given: length of a row of subway seats

Return: Mean and standard deviation of the fraction of occupied seats, assuming that people choose seats randomly but are unwilling to sit next to another person, until no acceptable seats remain.
"""


# Given a list of occupied seats, and number of total seats,
# return a list of available seats.
# getAvailable([0,2],6) -> [4,5]
def getAvailable(occupied,total):
    available = []
    for i in range(total):
        if i not in occupied and i-1 not in occupied and i+1 not in occupied:
            available.append(i)
    return available
def main():
    print getAvailable([0,2],6)

main()
