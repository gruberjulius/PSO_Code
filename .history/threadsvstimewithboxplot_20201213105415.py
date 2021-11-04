import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches


def plotJulius(Times, Threads, Deviations, name, xl, yl, i='', locationofstorage=''):
    """
    Plotting function for errorbars:
    Input are the arrays of Times and the Arrays of Threads, as well as the deviations
    Then the name of the function, the label for the x-Axes and for the y-Axes, 
    Optionally you can add the number of the plot as well as the saving location
    """
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    j = str(i)
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 17, 3, 4]) # scale the plotting
    plt.errorbar(Threads, Times, Deviations, linestyle='None', fmt='o', 
            color='black', ecolor='gray')
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    #set_xlim(1)
    #plt.legend( ['Medians of samples'], handles =Deviations)#, ['Standard Deviation'])
    plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    plt.title("plotting function number" + j)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) #xaxis.set_ticks([1.,2.,3.,10.])
    #plt.xaxis.set_ticks([1.,2.,3.,4., 5., 6., 7., 8., 9., 10., 11., 12.,, 13., 14., 15., 16.])
    #plt.ax1.plot(range(10))
    plt.savefig( locationofstorage + name + j + '.pdf')


def cleaning(Threads, Times, stepsize, numoffunc):
    """this functions "cleans" the data
    it taked in the big files with number of threads, the times, the standard step size used
    and the number of different functions one wants to plot
    it outputs a cleaned version, hence three arrays of size numoffunc
    the first is the numbering of threads, the second one with the averaged times for each threat step
    and the third one for the standarddeviations for each timestep
    """
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


def modifiedsortIndex(Index, Threads, Time, BestFitnessValue):
    """Function that will sort all the Elements according to their indixes
    Input will be the number of Threads, the times, the indces, the best fitness value
    Output will be the sorted arrays
    """
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


def modifiedsortIndexThreads(Index, Threads, Time, BestFitnessValue):
    """Function that will sort all the Elements according if their indices are equal to their thread number
    Input will be the number of Threads, the times, the indces, the best fitness value
    Output will be the sorted arrays according to Indixes and Arrays
    """
    n = len(Threads)
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
            if Index[j] == Index[j + 1] and Threads[j] > Threads[j + 1]:
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


def countDistinct(arr):
    """Takes an numpy array as an input 
    Ouput will be the number of different elements in it
    """  
    res = 1
    n = len(arr)
    # Pick all elements one by one 
    for i in range(1, n): 
        j = 0
        for j in range(i): 
            if (arr[i] == arr[j]): 
                break
  
        # If not printed earlier, then print it 
        if (i == j + 1): 
            res += 1      
    return res


def preperationwithdifferentindices(numberofdifferentfunctions, Index, Threads, Time, BestFitnessValue, stepsize):
    Numberofelements = len(Index)/ numberofdifferentfunctions
    Numberofelements = int(Numberofelements)
    i = 0
    for i in range(numberofdifferentfunctions):
        mThreads, mTimes, mDeviations = cleaning(Threads[i * Numberofelements: (i + 1) * Numberofelements], 
                Time[i * Numberofelements: (i + 1) * Numberofelements], stepsize, Numberofelements)
        plotJulius( mTimes, mThreads, mDeviations, "PlottingDifferentFunctions",
                "Number of Threads", "Times",i)


def main():
    """This is the critical function, the main function
    you will first need to input the location of the data set
    then you will need to add the number of runs, hence how many samples are done for each step
    the third input is the function name

    I ecpect that the data set has the following form:
    1. first index of the function
    2. the different outcomes for the time measuremtn
    3. The number of threads used in the current timestep
    4. the best fitness value of the function
    """
    location = input("The location of the file ")

    Index, Threads, Time, BestFitnessValue = np.loadtxt(
        #"results.dat", delimiter = ' ', unpack = True)
        location, delimiter=' ', unpack=True)

    # this will be 50 you guys said?
    stepsize = input("Number of runs for every thread ")
    stepsize = int(stepsize)
    # this will be 50 you guys said? kc: yes 50
    nameoffunc = input("What is the name of the function ")
    # stepsize = input("Number of runs for every thread") #this will be 50 you guys said?
    stepsize = 10  # int(stepsize)
    numberofthreads = int(len(Threads) / stepsize)
    # function that cleans the data
    numberofdifferentfunctions = countDistinct(Index)

    
    
    if(numberofdifferentfunctions == 1):
        mThreads, mTimes, mDeviations = cleaning(Threads, Time,
                                             stepsize, numberofthreads)
        plotJulius(mTimes, mThreads, mDeviations, nameoffunc,
               "Number of Threads", "Times")
    else:
        modifiedsortIndex(Index, Threads, Time, BestFitnessValue)
        modifiedsortIndexThreads(Index, Threads, Time, BestFitnessValue)
        preperationwithdifferentindices(numberofdifferentfunctions, Index, Threads, Time, BestFitnessValue, stepsize)


if __name__ == '__main__':
    main()
