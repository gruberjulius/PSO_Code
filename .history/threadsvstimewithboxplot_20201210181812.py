import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
#plt.style.use('seaborn-whitegrid')
# Times is the average Time needed
# Threads is the number of threads
# Deviation is the standard deviation
# name is the name of the plot
# xl stands for xlabel to label the x scale
# yl stands for ylabel to label the y scale


def plotJulius(Times, Threads, Deviations, name, xl, yl):
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    plt.axis([0, 17, 3, 4]) # scale the plotting
    plt.errorbar(Threads, Times, Deviations, linestyle='None', fmt='o', 
            color='black', ecolor='gray')
    plt.grid(linestyle='-',linewidth=1)
    #plt.legend( ['Medians of samples'], handles =Deviations)#, ['Standard Deviation'])
    plt.legend( handles = [Median_Plot, Deviations_Plot])
    plt.title(name)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0)

    plt.savefig(name + '.pdf')

# we need to give the function the number of samples per set
# numberoffuncs:number of threads input to the cleaning function i.e.15


# this functions "cleans" the data
# it taked in the big files with number of threads, the times, the standard step size used
# and the number of different functions one wants to plot
# it outputs a cleaned version, hence three arrays of size numoffunc
# the first is the numbering of threads, the second one with the averaged times for each threat step
# and the third one for the standarddeviations for each timestep
def cleaning(Threads, Times, stepsize, numoffunc):
    size = len(Threads)

    # get the number of different threads sizes used
    modifiedTimes = np.empty(numoffunc)
    modifiedDeviations = np.empty(numoffunc)
    modifiedThreads = np.arange(numoffunc)
    modifiedThreads = np.add(modifiedThreads, np.ones(numoffunc))

    i = 0

    while(i < numoffunc):
        modifiedTimes[i] = np.average(
            Times[(i * stepsize):((i * stepsize) + stepsize)])
        modifiedDeviations[i] = np.std(
            Times[(i * stepsize):((i * stepsize) + stepsize)])
        i += 1

    return modifiedThreads, modifiedTimes, modifiedDeviations


def modifiedsort(Index, Threads, Time, BestFitnessValue):
    n = len(Index)

    for i in range(n):
        # Create a flag that will allow the function to
        # terminate early if there's nothing left to sort
        already_sorted = True

        # Start looking at each item of the list one by one,
        # comparing it with its adjacent value. With each
        # iteration, the portion of the array that you look at
        # shrinks because the remaining items have already been
        # sorted.
        for j in range(n - i - 1):
            if Index[j] > Index[j + 1]:
                # If the item you're looking at is greater than its
                # adjacent value, then swap them
                # array[j], array[j + 1] = array[j + 1], array[j] #this seems to be the standard swapping function in python, should work
                Index[j], Index[j + 1] = Index[j + 1], Index[j]
                Threads[j], Threads[j + 1] = Threads[j + 1], Threads[j]
                Time[j], Time[j + 1] = Time[j + 1], Time[j]
                BestFitnessValue[j], BestFitnessValue[j +
                                                      1] = BestFitnessValue[j + 1], BestFitnessValue[j]
                # Since you had to swap two elements,
                # set the `already_sorted` flag to `False` so the
                # algorithm doesn't finish prematurely
                already_sorted = False
        # If there were no swaps during the last iteration,
        # the array is already sorted, and you can terminate
        if already_sorted:
            break

# This is the critical function, the main function
# you will first need to input the location of the data set
# then you will need to add the number of runs, hence how many samples are done for each step
# the third input is the function name

# I ecpect that the data set has the following form:
# 1. first index of the function
# 2. the different outcomes for the time measuremtn
# 3. The number of threads used in the current timestep
# 4. the best fitness value of the function


def main():
    location = input("The location of the file")

    Index, Threads, Time, BestFitnessValue = np.loadtxt(
        #"results.dat", delimiter = ' ', unpack = True)
        location, delimiter=' ', unpack=True)

    # this will be 50 you guys said?
    stepsize = input("Number of runs for every thread")
    stepsize = int(stepsize)
    # this will be 50 you guys said? kc: yes 50
    nameoffunc = input("What is the name of the function")
    # stepsize = input("Number of runs for every thread") #this will be 50 you guys said?
    stepsize = 10  # int(stepsize)
    numberofthreads = int(len(Threads) / stepsize)
    # function that cleans the data
    modifiedsort(Index, Threads, Time, BestFitnessValue)
    mThreads, mTimes, mDeviations = cleaning(Threads, Time,
                                             stepsize, numberofthreads)

    # plotting function
    plotJulius(mTimes, mThreads, mDeviations, nameoffunc,
               "Number of Threads", "Times")


if __name__ == '__main__':
    main()
