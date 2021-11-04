import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
import copy
import sys


def sumarrinv(arr):
    """
    function that calculates the sum of all the inverses of all the elements
    """
    size = len(arr)
    sum = 0 #introducing a summing variable
    i = 0 #introducing a counting variable
    while i < 0:
        if arr[i] > 0:
            sum += (1./arr[i])
            i += 1
        else:
            i += 1
    return sum


def hmean(arr):
    """calculates the harmonic mean of a functions
    the elements of the array have to be positive in order to get a accurate results
    """
    lenarr = len(arr) # get the length of the array
    sumarr = sumarrinv(arr) # get the sum of all the inverted elements
    if sumarr != 0 :
        return lenarr/sumarr # return the harmonic mean

#function that sorts the elements for the indicies
def modifiedsortIndex(Index, Threads, Time, BestFitnessValue):
    """
    Function that will sort all the Elements according to their indixes
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

#function to calculate the step size for the index and thread data
def numberofsamplesforIndexandThreads(Index, Thread):
    """
    function that returns the number of samples for each index and thread
    Input : Indexlist and Threadlist
    Output: Number of samples per Index and Threas
    """
    #Differentindicies = countDistinct(Index) should be the same
    Differentindicies = len(np.unique(Index))
    #DifferentThreads = countDistinct(Thread)
    DifferentThreads = len(np.unique(Thread))

    
    samplesize = len(Index)/( Differentindicies * DifferentThreads)
    #return samplesize
    return 50

#function that sorts elements in lexiographic order according to Index and Threads
def modifiedsortIndexThreads(Index, Threads, Time, BestFitnessValue):
    """
    Function that will sort all the Elements according if their indices are equal to their thread number
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

#bug fixed the function seperatly
def Sorting(Index, Thread, Time, BestFitnessValue):
    """
    function that calls the 2 sorting functions seperatly. Excluded this for readabilits
    """
    modifiedsortIndex(Index, Thread, Time, BestFitnessValue)
    print("keep your head up, we are currently sorting: ")
    modifiedsortIndexThreads(Index, Thread, Time, BestFitnessValue)

#currently on it
def ConvertDataDeviations(Time, PlottingDeviation, NumberofSamples):
    """
    Creating a seperate converting function for Deviations
    as they must be created seperatly
    Checked and it worked
    """
    i = 0 
    j = 0
    shape = np.shape(PlottingDeviation)
    a = PlottingDeviation.shape[0]
    b = PlottingDeviation.shape[1]
    while i < a:
        while j < b:
            begin = NumberofSamples*(i * b + j)
            begin = int(begin) 
            end = NumberofSamples*(i * b + j + 1)
            end = int(end)
            PlottingDeviation[i][j] = np.std(Time[begin: end])
            #the iteration is probably wrong, this is the wrong deviation
            j += 1
        j = 0
        i += 1
    return  PlottingDeviation


def ConvertDatahmean(Format, PlottingFormat, NumberofSamples):
    """
    Converts Data to different format to make it easier to plot
    """
    #need to hardcode this
    NumberofSamples = int(NumberofSamples)
    i = 0
    j = 0
    a = np.size(PlottingFormat, axis = 0)
    b = np.size(PlottingFormat, axis = 1)
    while i < a:
        while j < b:
            begin = NumberofSamples*(i * b + j)
            begin = int(begin) 
            end = NumberofSamples*(i * b + j + 1)
            end = int(end)
            if (len(np.unique(Format[begin:end])) == 1):
                PlottingFormat[i][j] = Format[begin]
            else: 
                PlottingFormat[i][j] = hmean(Format[begin: end])
            j += 1
        j = 0
        i += 1
    
    Format = PlottingFormat.copy()

#calculates the arithmetic mean
def ConvertDataamean(Format, PlottingFormat, NumberofSamples):
    """
    Converts Data to different format to make it easier to plot
    """
    #need to hardcode this
    NumberofSamples = int(NumberofSamples)
    i = 0
    j = 0
    a = np.size(PlottingFormat, axis = 0)
    b = np.size(PlottingFormat, axis = 1)
    while i < a:
        while j < b:
            begin = NumberofSamples*(i * a + j)
            begin = int(begin) 
            end = NumberofSamples*(i * a + j + 1)
            end = int(end)
            if (len(np.unique(Format[begin:end])) == 1):
                PlottingFormat[i][j] = Format[begin]
            else:
                PlottingFormat[i][j] = np.mean(Format[begin: end])
            #the iteration is probably wrong, this is the wrong deviation
            j += 1
        j = 0
        i += 1
    
    #Format = PlottingFormat.copy()

def ConvertData(Format, PlottingFormat):
    """
    Converts Data to different format to make it plotable in the end of
    """
    i = 0
    j = 0
    a = np.size(PlottingFormat, axis = 0)
    print(a)
    b = np.size(PlottingFormat, axis = 1)

    while i < a:
        while j < b:
            PlottingFormat[i][j] =  np.mean(Format[50 * (i  * b + j) : (i * b + j + 1)*50])
            j += 1
        j = 0
        i += 1

#need to fix the samplesize later on
#the non harmonic hence the arithmetic functions are called
#need to fix this later on
def Prepare(Index, Thread, Time, BestFitnessValue):
    print
    a = len(np.unique(Index))
    b = len(np.unique(Thread))

    #creating several 2 D arrays, which will be needed for the plotting
    PlottingIndex = np.empty((a, b)) # probaly convert this to integer
    PlottingThread =  np.empty((a, b)) # probaly convert this to integer
    PlottingTime =  np.empty((a, b))
    PlottingBestFitnessValue = np.empty((a, b))
    PlottingDeviation = np.empty((a, b))
    Samplesize = numberofsamplesforIndexandThreads(Index, Thread)    
    #the following equation must hold: NumberofSample * a * b
    #call different converting function to get data in the right format
    ConvertDataDeviations(Time, PlottingDeviation, Samplesize)
    ConvertData(Index, PlottingIndex)
    print(Index)
    print(PlottingIndex)
    ConvertData(Thread, PlottingThread)
    ConvertData(Time, PlottingTime)
    ConvertData(BestFitnessValue, PlottingBestFitnessValue)
    #we return all elements of the computation seperatly in order to call out plotting function later
    return PlottingIndex, PlottingThread, PlottingTime, PlottingBestFitnessValue, PlottingDeviation


def storeData(PlottingIndex, PlottingThread, PlottingDeviation, PlottingTime, PlottingBestFitnessValue , storagelocation):
    print(PlottingIndex, file=open(storagelocation, 'a'))
    print(PlottingThread, file=open(storagelocation, 'a'))
    print(PlottingTime , file=open(storagelocation, 'a'))
    print(PlottingDeviation, file=open(storagelocation, 'a'))
    print(PlottingBestFitnessValue, file=open(storagelocation, 'a'))
    #print(PlottingIndex)
    #print(PlottingThread)
    #print(PlottingTime)
    #print(PlottingBestFitnessValue )
    #print(PlottingDeviation)

def plot(Index, Thread, Time, Deviation, storagelocation ='', 
        titleoftheplot='', xl='Threads', yl = 'Times'):
    #I will create a seperate plot for every index and put them together in the end
    print("Currently plotting")
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 17, 3, 4]) # manually scale the plotting
    plt.title(titleoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) 
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = 0
    #print(len(np.unique(Index)))
    numberofdifferentindicies = len(np.unique(Index))
    ThreadDummy = Thread[i:(i+ 1)]
    TimeDummy = Time[i:(i + 1)]
    DeviationDummy = Deviation[i:(i + 1)]
    #print(ThreadDummy)
    #print(TimeDummy)
    #print(DeviationDummy)
    Time2Dummy = TimeDummy[0]
    Thread2Dummy = ThreadDummy[0]
    Deviation2Dummy = DeviationDummy[0]
    #print(Thread2Dummy)
    #print(Time2Dumm
    #print(Deviation2Dummy)
    #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy)
    #plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    #print("a")
    while i < 2:#numberofdifferentindicies:
        #print("b")
        ThreadDummy = Thread[i:(i+ 1)]
        DeviationDummy = Deviation[i:(i + 1)] #used this silly work arround
        TimeDummy = Time[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #print("c")
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        #print("d")
        i += 1
    plt.legend(loc='lower right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


def fillupp(Index, Thread, Time, Bestfitnessvalue):
    print("Currently filling up")
    numberofdifferentindicies = len(np.unique(Index)) # 6
    print(numberofdifferentindicies)
    numberofdifferentthreads = len(np.unique(Thread)) # 13 
    print(numberofdifferentthreads)
    numberofcurrentdata = 0 #np.zero((numberofdifferentindicies, numberofdifferentthreads))
    expectedsize = numberofdifferentindicies * numberofdifferentthreads * 50 # size the array should have
    newIndex = np.empty(expectedsize)
    newThread = np.empty(expectedsize)
    newTime = np.empty(expectedsize)
    newBestFV = np.empty(expectedsize)
    print(len(newIndex))
    print(len(Index))

    i = 0
    j = 0
    k1 = 0
    k2 = 0
    while k1 < len(Index) and k2 < expectedsize:
        while i < numberofdifferentindicies:
            numberofcurrentdata = 0
            while j < numberofdifferentthreads:

                if Index[k1] == i and Thread[k1] == j:
                    numberofcurrentdata +=1
                    if numberofcurrentdata <= 50:                        
                        newIndex[k1] = Index[k1]
                        newTime[k1] = Time[k1]
                        newThread[k1] = Thread[k1]
                        newBestFV[k1] = Bestfitnessvalue[k1]
                j += 1
            if numberofcurrentdata < 50:
                l = 50 - numberofcurrentdata 
                k2 -= 1
                for m in range(0,l):
                    k2 += 1
                    newIndex[k2] = Index[k1 -m]
                    newTime[k2] = Time[k1- m]
                    newThread[k2] = Thread[k1- m]
                    newBestFV[k2] = Bestfitnessvalue[k1- m]
                    #print(k2)
            
                #print(Index[i])
                #print(Thread[j])
            j = 0
            i += 1
        i = 0
        k2 += 1
        k1 += 1     
    #print(k2)  
    if len(newTime) == len(Time):
        print("you arre lucky bastards!!!!!")
    return newIndex, newThread, newTime, newBestFV      


def preclean(nameoffile):
    print("Currently precleaning")
    with open(nameoffile, "r") as f:
        lines = f.readlines()
    with open(nameoffile, "w") as f:
        for line in lines:
            if len(line.strip("\n").split(" ")) == 4:
                f.write(line) 


def main():
    """
    I radically changed the build up of my function. I propose the new following build up.
    1. A seperate readin function for readability
    2. A Sorting function that call the 2 different sorting patterns at the same time.
    3. A preperation function that resorts the storage of the function so that it will be easier to access for plotting 
    """
    nameoffile = input("What is the location of your file! ")
    #preclean(nameoffile)
    oIndex , oThread , oTime , oBestFitnessValue = np.loadtxt( nameoffile , delimiter=' ', unpack=True)
    Sorting(oIndex, oThread, oTime, oBestFitnessValue)
    Index, Thread, Time, BestFitnessValue =fillupp(oIndex, oThread, oTime, oBestFitnessValue)
    PlottingIndex, PlottingThread, PlottingTime, PlottingBestFitnessValue, PlottingDeviation =\
                Deviations = Prepare(Index, Thread, Time, BestFitnessValue) 

    #storagelocation = input("Input the storage location, where you want to store your files? ")
    #nameoftheplot = input("What is the name of the plot")
    #plot(PlottingIndex, PlottingThread, PlottingTime, PlottingDeviation)#, storagelocation, nameoftheplot)
    #location and name of the plot will be added later on
    #storagelocation = "example.npy"
    #storeData(PlottingIndex, PlottingThread, PlottingTime,  PlottingDeviation,  
    #        PlottingBestFitnessValue, storagelocation)
    titleofplot = input("What is the title of the plot? ")
    plot(PlottingIndex, PlottingThread, PlottingTime, PlottingDeviation, " ", titleofplot)
    print(PlottingIndex)
    print(PlottingThread)
    print(PlottingTime)
    print(PlottingBestFitnessValue )
    print(PlottingDeviation)


if __name__ == '__main__':
    main()
