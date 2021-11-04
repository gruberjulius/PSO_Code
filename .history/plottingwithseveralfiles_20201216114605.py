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
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = 0
    #print(len(np.unique(Index)))
    numberofdifferentindicies = len(np.unique(Index))
    #ThreadDummy = Thread[i:(i+ 1)]
    #TimeDummy = Time[i:(i + 1)]
    #DeviationDummy = Deviation[i:(i + 1)]
    #print(ThreadDummy)
    #print(TimeDummy)
    #print(DeviationDummy)
    #Time2Dummy = TimeDummy[0]
    #Thread2Dummy = ThreadDummy[0]
    #Deviation2Dummy = DeviationDummy[0]
    #print(Thread2Dummy)
    #print(Time2Dumm
    #print(Deviation2Dummy)
    #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy)
    #plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    #print("a")
    while i < numberofdifferentindicies:
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
    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


def plotwith2inputfiles(Index1, Index2, Thread1, Thread2, Time1, Time2, Deviation1, Deviation2, storagelocation ='', 
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
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    
    i = 0
    #print(len(np.unique(Index)))
    numberofdifferentindicies = len(np.unique(Index))
    #ThreadDummy = Thread[i:(i+ 1)]
    #TimeDummy = Time[i:(i + 1)]
    #DeviationDummy = Deviation[i:(i + 1)]
    #print(ThreadDummy)
    #print(TimeDummy)
    #print(DeviationDummy)
    #Time2Dummy = TimeDummy[0]
    #Thread2Dummy = ThreadDummy[0]
    #Deviation2Dummy = DeviationDummy[0]
    #print(Thread2Dummy)
    #print(Time2Dumm
    #print(Deviation2Dummy)
    #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy)
    #plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    #print("a")
    while i < numberofdifferentindicies:
        #print("b")
        ThreadDummy = Thread1[i:(i+ 1)]
        DeviationDummy = Deviation1[i:(i + 1)] #used this silly work arround
        TimeDummy = Time1[i:(i + 1)]
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
    i = 0
    while i < numberofdifferentindicies :
        #print("b")
        ThreadDummy = Thread2[i:(i+ 1)]
        DeviationDummy = Deviation2[i:(i + 1)] #used this silly work arround
        TimeDummy = Time2[i:(i + 1)]
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

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


def plotwith3inputfiles(Index1, Index2,Index3,  Thread1, Thread2, Thread3,  Time1, Time2, Time3, Deviation1, Deviation2, Deviation3, 
        storagelocation ='', titleoftheplot='', xl='Threads', yl = 'Times'):
    #I will create a seperate plot for every index and put them together in the end
    print("Currently plotting")
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
    #print(len(np.unique(Index)))
    numberofdifferentindicies = len(np.unique(Index))
    #ThreadDummy = Thread[i:(i+ 1)]
    #TimeDummy = Time[i:(i + 1)]
    #DeviationDummy = Deviation[i:(i + 1)]
    #print(ThreadDummy)
    #print(TimeDummy)
    #print(DeviationDummy)
    #Time2Dummy = TimeDummy[0]
    #Thread2Dummy = ThreadDummy[0]
    #Deviation2Dummy = DeviationDummy[0]
    #print(Thread2Dummy)
    #print(Time2Dumm
    #print(Deviation2Dummy)
    #plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy)
    #plt.legend( handles = [Median_Plot, Deviations_Plot]) # can add more names to the legend now
    #print("a")
    while i < numberofdifferentindicies:
        #print("b")
        ThreadDummy = Thread1[i:(i+ 1)]
        DeviationDummy = Deviation1[i:(i + 1)] #used this silly work arround
        TimeDummy = Time1[i:(i + 1)]
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

    i = 0 #setting the index variable back to zero
    while i < numberofdifferentindicies :
        #print("b")
        ThreadDummy = Thread2[i:(i+ 1)]
        DeviationDummy = Deviation2[i:(i + 1)] #used this silly work arround
        TimeDummy = Time2[i:(i + 1)]
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
    
    i = 0
    while i < numberofdifferentindicies :
        #print("b")
        ThreadDummy = Thread3[i:(i+ 1)]
        DeviationDummy = Deviation3[i:(i + 1)] #used this silly work arround
        TimeDummy = Time3[i:(i + 1)]
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

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + titleoftheplot +'.pdf' )


def fillupp(Index, Thread, Time, Bestfitnessvalue):
    print("Currently filling up")
    numberofdifferentindicies = len(np.unique(Index)) # 6
    print(numberofdifferentindicies)
    numberofdifferentthreads = len(np.unique(Thread)) # 13 
    print(numberofdifferentthreads)
    numberofcurrentdata = 0 #np.zero((numberofdifferentindicies, numberofdifferentthreads))
    expectedsize = numberofdifferentindicies * numberofdifferentthreads * 50 # size the array should have
    if(expectedsize == len(Index)):
        return 0
    newIndex = np.zeros(expectedsize)
    newThread = np.zeros(expectedsize)
    newTime = np.zeros(expectedsize)
    newBestFV = np.zeros(expectedsize)
    print(newIndex)
    numberofcurrentdata = np.zeros(numberofdifferentindicies * numberofdifferentthreads)
    print(len(newIndex))
    print(len(Index))
    ThreadId = np.unique(Thread)
    print(ThreadId)
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
    print(numberofcurrentdata)
    print(np.sum(numberofcurrentdata))
    #print(k2)  
    i = 0
    j = 0
    k1 = 0
    k2 = 0
    kold = 0
    while k1 < len(Index) and k2 < expectedsize:
        while i < numberofdifferentindicies:
            while j < numberofdifferentthreads:
                currentnumber = int(numberofcurrentdata[i* numberofdifferentthreads + j])
                print(currentnumber)
                #if Index[k1] == i and Thread[k1] == j:
                #    numberofcurrentdata[i* numberofdifferentindicies + j] +=1 #hopefully this is correctl working
                if currentnumber == 50: 
                    for o in range(0, 50):  
                        #print(k1)
                        #print(Index[k1])                     
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
                        newIndex[k2] = Index[k1 - m]
                        newTime[k2] = Time[k1 - m]
                        newThread[k2] = Thread[k1 - m]
                        newBestFV[k2] = Bestfitnessvalue[k1 - m]
                        k2 += 1
                j += 1
            j = 0
            i += 1
        i = 0
    print(numberofcurrentdata)
    print("new index")
    print(newIndex)
    print("new Thread")
    print( newThread)
    print( newTime)
    print( newBestFV)    
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
        plotwith2inputfiles(PlottingIndex1, PlottingIndex2, PlottingThread1, PlottingThread2, PlottingTime1, PlottingTime2, 
            PlottingDeviation1, PlottingDeviation2, titleofplot, titleofplot)
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
        plot(PlottingIndex, PlottingThread, PlottingTime, PlottingDeviation, " ", titleofplot)



if __name__ == '__main__':
    main()
