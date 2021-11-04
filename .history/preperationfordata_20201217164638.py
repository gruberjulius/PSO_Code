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
    #print("keep your head up, we are currently sorting: ")
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
            PlottingDeviation[i][j] = hmean(Time[begin: end])
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



def fillupp(Index, Thread, Time, Bestfitnessvalue):
    #print("Currently filling up")
    numberofdifferentindicies = len(np.unique(Index)) # 6
    numberofdifferentthreads = len(np.unique(Thread)) # 13 
    numberofcurrentdata = 0 #np.zero((numberofdifferentindicies, numberofdifferentthreads))
    expectedsize = numberofdifferentindicies * numberofdifferentthreads * 50 # size the array should have
    #print(expectedsize)
    if(expectedsize == len(Index)):
        return 0
    newIndex = np.zeros(expectedsize)
    newThread = np.zeros(expectedsize)
    newTime = np.zeros(expectedsize)
    newBestFV = np.zeros(expectedsize)
    numberofcurrentdata = np.zeros(numberofdifferentindicies * numberofdifferentthreads)
    ThreadId = np.unique(Thread)
    i = 0
    j = 0
    k1 = 0
    k2 = 0
    kold = 0
    while k1 < len(Index): #and k2 < expectedsize:
        while i < numberofdifferentindicies:
            while j < numberofdifferentthreads:
                if Index[k1] == i and Thread[k1] == ThreadId[j]:
                    numberofcurrentdata[i* numberofdifferentthreads + j] +=1 #hopefully this is correctl working  
                j +=1
            j = 0
            i += 1
        i = 0
        #k2 += 1
        k1 += 1     
    i = 0
    j = 0
    k1 = 0
    k2 = 0
    kold = 0
    #print(numberofcurrentdata)
    while k1 < len(Index) and k2 < expectedsize:
        while i < numberofdifferentindicies:
            while j < numberofdifferentthreads:
                currentnumber = int(numberofcurrentdata[i* numberofdifferentthreads + j])
                #if Index[k1] == i and Thread[k1] == j:
                #    numberofcurrentdata[i* numberofdifferentindicies + j] +=1 #hopefully this is correctl working
                if currentnumber == 50: 
                    for o in range(0, 50):  
                        newIndex[k2] = Index[k1]
                        newTime[k2] = Time[k1]
                        newThread[k2] = Thread[k1]
                        newBestFV[k2] = Bestfitnessvalue[k1]
                        k1 += 1
                        k2 += 1 
                elif currentnumber > 50:
                    for q in range(0, currentnumber):
                        if q < 50:
                            newIndex[k2] = Index[k1]
                            newTime[k2] = Time[k1]
                            newThread[k2] = Thread[k1]
                            newBestFV[k2] = Bestfitnessvalue[k1]
                            k1 += 1
                            k2 += 1
                        else:
                            k1 += 1
                else : # currentnumber < 50
                    l = 50 - currentnumber                         
                    for q in range( 0, currentnumber):
                        newIndex[k2] = Index[k1]
                        newTime[k2] = Time[k1]
                        newThread[k2] = Thread[k1]
                        newBestFV[k2] = Bestfitnessvalue[k1]
                        k1 += 1
                        k2 += 1
                    for m in range(0,l):
                        newIndex[k2] = Index[k1 - 1 ]
                        newTime[k2] = Time[k1 - 1 ]
                        newThread[k2] = Thread[k1 - 1  ]
                        newBestFV[k2] = Bestfitnessvalue[k1 -1 ]
                        k2 += 1
                        
                    
                j += 1
            j = 0
            i += 1
        i = 0
    return newIndex, newThread, newTime, newBestFV      


def preclean(nameoffile):
    #print("Currently precleaning")
    with open(nameoffile, "r") as f:
        lines = f.readlines()
    with open(nameoffile, "w") as f:
        for line in lines:
            if len(line.strip("\n").split(" ")) == 4:
                f.write(line) 


def SaveData(Index, Thread, Time, BestFitnessValue, nameoffile):
    numberofsamples = len(Index)
    file = open(nameoffile, "w")
    i = 0 #counting variable
    while i < numberofsamples:
        line = str(Index[i]) + " " + str(Thread[i]) + " " +  str(Time[i]) + " " +  str(BestFitnessValue[i]) + "\n"
        file.write(line)
        i += 1

    file.close()


def main():
    """
    I radically changed the build up of my function. I propose the new following build up.
    1. A seperate readin function for readability
    2. A Sorting function that call the 2 different sorting patterns at the same time.
    3. A preperation function that resorts the storage of the function so that it will be easier to access for plotting 
    """
    #nameoffile = input("What is the location of your file! ")
    nameoffile = input("What is the name of your file? ")
    preclean(nameoffile)
    oIndex,oThread, oTime , oBestFitnessValue = np.loadtxt( nameoffile , delimiter=' ', unpack=True)
    Sorting(oIndex, oThread, oTime, oBestFitnessValue)
    Index, Thread, Time, BestFitnessValue = fillupp(oIndex, oThread, oTime, oBestFitnessValue)
    #numberofsamples = len(Index)
    #file = open(nameoffile, "w")
    #i = 0 #counting variable
    #while i < numberofsamples:
    #    line = str(Index[i]) + " " + str(Thread[i]) + " " +  str(Time[i]) + " " +  str(BestFitnessValue[i]) + "\n"
    #    file.write(line)
    #    i += 1

    #file.close()
    SaveData(Index, Thread, Time, BestFitnessValue, nameoffile)

if __name__ == '__main__':
    main()