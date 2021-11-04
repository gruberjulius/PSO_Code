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


def plot(Index, Thread, Time, Deviation, storagelocation ='', 
        titleoftheplot='', xl='Threads', yl = 'Times'):
    #I will create a seperate plot for every index and put them together in the end
    #print("Currently plotting")
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 17, 3, 4]) # manually scale the plotting
    plt.title(titleoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) 
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = 0
    numberofdifferentindicies = len(np.unique(Index))

    while i < numberofdifferentindicies:
        ThreadDummy = Thread[i:(i+ 1)]
        DeviationDummy = Deviation[i:(i + 1)] #used this silly work arround
        TimeDummy = Time[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1
    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


def plotwith2inputfiles( Index1, Index2, Thread1, Thread2, Time1, Time2, Deviation1, Deviation2, storagelocation , 
        titleoftheplot, startindex , endindex,  xl='Threads', yl = 'Times'):
    #I will create a seperate plot for every index and put them together in the end
    #print("Currently plotting")
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 17, 3, 4]) # manually scale the plotting
    plt.title(titleoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) 
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = int(startindex)
    k = int(endindex)
    #plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    while i < k :#numberofdifferentindicies1:
        ThreadDummy = Thread1[i:(i+ 1)]
        DeviationDummy = Deviation1[i:(i + 1)] #used this silly work arround
        TimeDummy = Time1[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o',  label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1
    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies2 :
        ThreadDummy = Thread2[i:(i+ 1)]
        DeviationDummy = Deviation2[i:(i + 1)] #used this silly work arround
        TimeDummy = Time2[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='x', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


def plotwith3inputfiles(Index1, Index2,Index3,  Thread1, Thread2, Thread3,  Time1, Time2, Time3, Deviation1, Deviation2, Deviation3, 
        storagelocation ='', titleoftheplot='', xl='Threads', yl = 'Times'):
    #I will create a seperate plot for every index and put them together in the end
    #("Currently plotting")
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 17, 3, 4]) # manually scale the plotting
    plt.title(titleoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) 
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = 0
    numberofdifferentindicies = len(np.unique(Index))
    while i < numberofdifferentindicies:
        ThreadDummy = Thread1[i:(i+ 1)]
        DeviationDummy = Deviation1[i:(i + 1)] #used this silly work arround
        TimeDummy = Time1[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1

    i = 0 #setting the index variable back to zero
    while i < numberofdifferentindicies :
        ThreadDummy = Thread2[i:(i+ 1)]
        DeviationDummy = Deviation2[i:(i + 1)] #used this silly work arround
        TimeDummy = Time2[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1
    
    i = 0
    while i < numberofdifferentindicies :
        ThreadDummy = Thread3[i:(i+ 1)]
        DeviationDummy = Deviation3[i:(i + 1)] #used this silly work arround
        TimeDummy = Time3[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )

#you need to manually adjust the legend
def plotwith7inputfiles(Index1, Index2,Index3,Index4, Index5, Index6, Index7, Thread1, Thread2, Thread3, Thread4,
        Thread5, Thread6, Thread7, Time1, Time2, Time3, Time4, Time5, Time6, Time7, Deviation1, Deviation2, 
        Deviation3, Deviation4, Deviation5, Deviation6, Deviation7, 
        storagelocation , titleoftheplot, startindex, endindex, xl='Threads', yl = 'Times'):
    #I will create a seperate plot for every index and put them together in the end
    #print("Currently plotting")
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 50, 0, 50]) # manually scale the plotting
    plt.title(titleoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) 
    #Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    #Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = int(startindex)
    k = int(endindex)
    #plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    while i < k :#numberofdifferentindicies1:
        ThreadDummy = Thread1[i:(i+ 1)]
        DeviationDummy = Deviation1[i:(i + 1)] #used this silly work arround
        TimeDummy = Time1[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='o',  label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        #plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        plt.plot(Thread2Dummy, Time2Dummy)
        i += 1

    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies2 :
        ThreadDummy = Thread2[i:(i+ 1)]
        DeviationDummy = Deviation2[i:(i + 1)] #used this silly work arround
        TimeDummy = Time2[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='x', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1
    
    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies3 :
        ThreadDummy = Thread3[i:(i+ 1)]
        DeviationDummy = Deviation3[i:(i + 1)] #used this silly work arround
        TimeDummy = Time3[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='^', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1

    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies4 :
        ThreadDummy = Thread4[i:(i+ 1)]
        DeviationDummy = Deviation4[i:(i + 1)] #used this silly work arround
        TimeDummy = Time4[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='H', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1    
    
    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies5 :
        ThreadDummy = Thread5[i:(i+ 1)]
        DeviationDummy = Deviation5[i:(i + 1)] #used this silly work arround
        TimeDummy = Time5[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='s', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1
    
    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies6 :
        ThreadDummy = Thread6[i:(i+ 1)]
        DeviationDummy = Deviation6[i:(i + 1)] #used this silly work arround
        TimeDummy = Time6[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='D', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1
    
    i = int(startindex)
    k = int(endindex)
    while i < k :#numberofdifferentindicies6 :
        ThreadDummy = Thread7[i:(i+ 1)]
        DeviationDummy = Deviation7[i:(i + 1)] #used this silly work arround
        TimeDummy = Time7[i:(i + 1)]
        Time2Dummy = TimeDummy[0]
        Thread2Dummy = ThreadDummy[0]
        Deviation2Dummy = DeviationDummy[0]
        j = str(i)
        #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, linestyle='None', fmt='2', label ='This is the plot of function ' + j)#, 
            #color='black', ecolor='gray')   
        plt.plot(Thread2Dummy, Time2Dummy)#, label = 'Function with index ' + j)
        i += 1  

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


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


def main():
    """
    I radically changed the build up of my function. I propose the new following build up.
    1. A seperate readin function for readability
    2. A Sorting function that call the 2 different sorting patterns at the same time.
    3. A preperation function that resorts the storage of the function so that it will be easier to access for plotting 
    """
    #nameoffile = input("What is the location of your file! ")
    placetogo = int(input("How many files do you want to plot together? "))
    if placetogo == 1:
        nameoffile = input("What is the name of your file? ")
        preclean(nameoffile)
        oIndex,oThread, oTime , oBestFitnessValue = np.loadtxt( nameoffile , delimiter=' ', unpack=True)
        Sorting(oIndex, oThread, oTime, oBestFitnessValue)
        Index, Thread, Time, BestFitnessValue = fillupp(oIndex, oThread, oTime, oBestFitnessValue)
        PlottingIndex, PlottingThread, PlottingTime, PlottingBestFitnessValue, PlottingDeviation =\
                    Deviations = Prepare(Index, Thread, Time, BestFitnessValue)
        titleofplot = input("What is the title of the plot? ")
        plot(PlottingIndex, PlottingThread, PlottingTime, PlottingDeviation, " ", titleofplot) 
    elif placetogo == 2:
        nameoffile1 = input("What is the name of your file 1? ")
        nameoffile2 = input("What is the name of your file 2? ")
        preclean(nameoffile1)
        preclean(nameoffile2)
        oIndex1, oThread1, oTime1 , oBestFitnessValue1 = np.loadtxt( nameoffile1 , delimiter=' ', unpack=True)
        oIndex2, oThread2, oTime2 , oBestFitnessValue2 = np.loadtxt( nameoffile2 , delimiter=' ', unpack=True)
        Sorting(oIndex1, oThread1, oTime1, oBestFitnessValue1)
        Sorting(oIndex2, oThread2, oTime2, oBestFitnessValue2)
        Index1, Thread1, Time1, BestFitnessValue1 = fillupp(oIndex1, oThread1, oTime1, oBestFitnessValue1)
        Index2, Thread2, Time2, BestFitnessValue2 = fillupp(oIndex2, oThread2, oTime2, oBestFitnessValue2)
        PlottingIndex1, PlottingThread1, PlottingTime1, PlottingBestFitnessValue1, PlottingDeviation1 =\
                    Deviations = Prepare(Index1, Thread1, Time1, BestFitnessValue1) 
        PlottingIndex2, PlottingThread2, PlottingTime2, PlottingBestFitnessValue2, PlottingDeviation2 =\
                    Deviations = Prepare(Index2, Thread2, Time2, BestFitnessValue2)
        titleofplot = input("What is the title of the plot? ") 
        startindex = int(input("What is the starting index of the functions you want to plot: "))
        endindex = int(input("what is the end index of functions U use? ")) # exklusiv
        xlabel = input("how do you want to label your x axes: ")
        ylabel = input("how do you want to label your y axes: ")
        plotwith2inputfiles(PlottingIndex1, PlottingIndex2, PlottingThread1, PlottingThread2, PlottingTime1, PlottingTime2, 
            PlottingDeviation1, PlottingDeviation2, titleofplot, titleofplot, startindex, endindex, xlabel, ylabel )
    elif placetogo == 3:
        nameoffile1 = input("What is the name of your file 1? ")
        nameoffile2 = input("What is the name of your file 2? ")
        nameoffile3 = input("What is the name of your file 3? ")
        preclean(nameoffile1)
        preclean(nameoffile2)
        preclean(nameoffile3)

        oIndex1, oThread1, oTime1 , oBestFitnessValue1 = np.loadtxt( nameoffile1 , delimiter=' ', unpack=True)
        oIndex2, oThread2, oTime2 , oBestFitnessValue2 = np.loadtxt( nameoffile2 , delimiter=' ', unpack=True)
        oIndex3, oThread3, oTime3 , oBestFitnessValue3 = np.loadtxt( nameoffile3 , delimiter=' ', unpack=True)
        Sorting(oIndex1, oThread1, oTime1, oBestFitnessValue1)
        Sorting(oIndex2, oThread2, oTime2, oBestFitnessValue2)
        Sorting(oIndex3, oThread3, oTime3, oBestFitnessValue3)
        Index1, Thread1, Time1, BestFitnessValue1 = fillupp(oIndex1, oThread1, oTime1, oBestFitnessValue1)
        Index2, Thread2, Time2, BestFitnessValue2 = fillupp(oIndex2, oThread2, oTime2, oBestFitnessValue2)
        Index3, Thread3, Time3, BestFitnessValue3 = fillupp(oIndex3, oThread3, oTime3, oBestFitnessValue3)
        PlottingIndex1, PlottingThread1, PlottingTime1, PlottingBestFitnessValue1, PlottingDeviation1 =\
                    Deviations = Prepare(Index1, Thread1, Time1, BestFitnessValue1) 
        PlottingIndex2, PlottingThread2, PlottingTime2, PlottingBestFitnessValue2, PlottingDeviation2 =\
                    Deviations = Prepare(Index2, Thread2, Time2, BestFitnessValue2)
        PlottingIndex3, PlottingThread3, PlottingTime3, PlottingBestFitnessValue3, PlottingDeviation3 =\
                    Deviations = Prepare(Index3, Thread3, Time3, BestFitnessValue3)
        titleofplot = input("What is the title of the plot? ") 
        plotwith3inputfiles(PlottingIndex1, PlottingIndex2, PlottingIndex3, PlottingThread1, 
            PlottingThread2, PlottingThread3, PlottingTime1, PlottingTime2, PlottingTime3, PlottingDeviation1, 
            PlottingDeviation2,PlottingDeviation3,  titleoftheplot, titleoftheplot)


            #storagelocation = input("Input the storage location, where you want to store your files? ")
            #nameoftheplot = input("What is the name of the plot")
            #plot(PlottingIndex, PlottingThread, PlottingTime, PlottingDeviation)#, storagelocation, nameoftheplot)
            #location and name of the plot will be added later on
            #storagelocation = "example.npy"
            #storeData(PlottingIndex, PlottingThread, PlottingTime,  PlottingDeviation,  
            #        PlottingBestFitnessValue, storagelocation)
        titleofplot = input("What is the title of the plot? ")
        plotwith3inputfiles(Index1, Index2, Index3, Thread1, Thread2, Thread3, Time1, Time2, Time3, Deviation1, Deviation2, Deviation3, storagelocation, titleoftheplot) (PlottingIndex, PlottingThread, PlottingTime, PlottingDeviation, " ", titleofplot)
    elif placetogo == 7:
        nameoffile1 = input("What is the name of your file 1? ")
        nameoffile2 = input("What is the name of your file 2? ")
        nameoffile3 = input("What is the name of your file 3? ")
        nameoffile4 = input("What is the name of your file 4? ")
        nameoffile5 = input("What is the name of your file 5? ")
        nameoffile6 = input("What is the name of your file 6? ")
        nameoffile7 = input("What is the name of your file 7? ")
        preclean(nameoffile1)
        preclean(nameoffile2)
        preclean(nameoffile3)
        preclean(nameoffile4)
        preclean(nameoffile5)
        preclean(nameoffile6)
        preclean(nameoffile7)
        oIndex1, oThread1, oTime1 , oBestFitnessValue1 = np.loadtxt( nameoffile1 , delimiter=' ', unpack=True)
        oIndex2, oThread2, oTime2 , oBestFitnessValue2 = np.loadtxt( nameoffile2 , delimiter=' ', unpack=True)
        oIndex3, oThread3, oTime3 , oBestFitnessValue3 = np.loadtxt( nameoffile3 , delimiter=' ', unpack=True)
        oIndex4, oThread4, oTime4 , oBestFitnessValue4 = np.loadtxt( nameoffile4 , delimiter=' ', unpack=True)
        oIndex5, oThread5, oTime5 , oBestFitnessValue5 = np.loadtxt( nameoffile5 , delimiter=' ', unpack=True)
        oIndex6, oThread6, oTime6 , oBestFitnessValue6 = np.loadtxt( nameoffile6 , delimiter=' ', unpack=True)
        oIndex7, oThread7, oTime7 , oBestFitnessValue7 = np.loadtxt( nameoffile7 , delimiter=' ', unpack=True) 
        Sorting(oIndex1, oThread1, oTime1, oBestFitnessValue1)
        Sorting(oIndex2, oThread2, oTime2, oBestFitnessValue2)
        Sorting(oIndex3, oThread3, oTime3, oBestFitnessValue3)
        Sorting(oIndex4, oThread4, oTime4, oBestFitnessValue4)
        Sorting(oIndex5, oThread5, oTime5, oBestFitnessValue5)
        Sorting(oIndex6, oThread6, oTime6, oBestFitnessValue6)
        Sorting(oIndex7, oThread7, oTime7, oBestFitnessValue7)
        Index1, Thread1, Time1, BestFitnessValue1 = fillupp(oIndex1, oThread1, oTime1, oBestFitnessValue1)
        Index2, Thread2, Time2, BestFitnessValue2 = fillupp(oIndex2, oThread2, oTime2, oBestFitnessValue2)
        Index3, Thread3, Time3, BestFitnessValue3 = fillupp(oIndex3, oThread3, oTime3, oBestFitnessValue3)
        Index4, Thread4, Time4, BestFitnessValue4 = fillupp(oIndex4, oThread4, oTime4, oBestFitnessValue4)
        Index5, Thread5, Time5, BestFitnessValue5 = fillupp(oIndex5, oThread5, oTime5, oBestFitnessValue5)
        Index6, Thread6, Time6, BestFitnessValue6 = fillupp(oIndex6, oThread6, oTime6, oBestFitnessValue6)
        Index7, Thread7, Time7, BestFitnessValue7 = fillupp(oIndex7, oThread7, oTime7, oBestFitnessValue7)
        PlottingIndex1, PlottingThread1, PlottingTime1, PlottingBestFitnessValue1, PlottingDeviation1 =\
                    Deviations = Prepare(Index1, Thread1, Time1, BestFitnessValue1) 
        PlottingIndex2, PlottingThread2, PlottingTime2, PlottingBestFitnessValue2, PlottingDeviation2 =\
                    Deviations = Prepare(Index2, Thread2, Time2, BestFitnessValue2)
        PlottingIndex3, PlottingThread3, PlottingTime3, PlottingBestFitnessValue3, PlottingDeviation3 =\
                    Deviations = Prepare(Index3, Thread3, Time3, BestFitnessValue3)
        PlottingIndex4, PlottingThread4, PlottingTime4, PlottingBestFitnessValue4, PlottingDeviation4 =\
                    Deviations = Prepare(Index4, Thread4, Time4, BestFitnessValue4)
        PlottingIndex5, PlottingThread5, PlottingTime5, PlottingBestFitnessValue5, PlottingDeviation5 =\
                    Deviations = Prepare(Index5, Thread5, Time5, BestFitnessValue5)
        PlottingIndex6, PlottingThread6, PlottingTime6, PlottingBestFitnessValue6, PlottingDeviation6 =\
                    Deviations = Prepare(Index6, Thread6, Time6, BestFitnessValue6)
        PlottingIndex7, PlottingThread7, PlottingTime7, PlottingBestFitnessValue7, PlottingDeviation7 =\
                    Deviations = Prepare(Index7, Thread7, Time7, BestFitnessValue7)


        titleofplot = input("What is the title of the plot? ") 
        startindex = int(input("What is the starting index of the functions you want to plot: "))
        endindex = int(input("what is the end index of functions U use? ")) # exklusiv
        xlabel = input("how do you want to label your x axes: ")
        ylabel = input("how do you want to label your y axes: ")
        storagelocation = input("What is storagelocation")
        plotwith7inputfiles( PlottingIndex1, PlottingIndex2, PlottingIndex3, PlottingIndex4, PlottingIndex5, PlottingIndex6, PlottingIndex7, 
            PlottingThread1, PlottingThread2, PlottingThread3, PlottingThread4, PlottingThread5, PlottingThread6, PlottingThread7, 
            PlottingTime1, PlottingTime2, PlottingTime3, PlottingTime4, PlottingTime5, PlottingTime6, PlottingTime7, PlottingBestFitnessValue1, 
            PlottingBestFitnessValue2, PlottingBestFitnessValue3, PlottingBestFitnessValue4, PlottingBestFitnessValue5, PlottingBestFitnessValue6, 
            PlottingBestFitnessValue7, storagelocation, titleofplot, startindex, endindex, xlabel, ylabel)

#def plotwith7inputfiles(Index1, Index2,Index3,Index4, Index5, Index6, Index7, Thread1, Thread2, Thread3, Thread4,
#        Thread5, Thread6, Thread7, Time1, Time2, Time3, Time4, Time5, Time6, Time7, Deviation1, Deviation2, 
#        Deviation3, Deviation4, Deviation5, Deviation6, Deviation7, 
#        storagelocation , titleoftheplot, startindex, endindex, xl='Threads', yl = 'Times'):

if __name__ == '__main__':
    main()
