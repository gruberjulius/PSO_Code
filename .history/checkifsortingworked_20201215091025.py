import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches

def hmean(arr):
    """calculates the harmonic mean of a functions
    the elements of the array have to be positive in order to get a accurate results
    """
    inverse = 1./arr #invert the array elements
    lenarr = len(inverse) # get the length of the array
    sumarr = sum(inverse) # get the sum of all the inverted elements
    if sumarr != 0 :
        return lenarr/sumarr # return the harmonic mean


def countDistinct(arr):

    """
    Takes an numpy array as an input 
    Ouput will be the number of different elements and the different elements in a numpy array
    Tested it works fine
    """  
    res = 0
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
    return samplesize

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
    modifiedsortIndexThreads(Index, Thread, Time, BestFitnessValue)

#currently on it
def ConvertDataDeviations(Time, PlottingDeviation, NumberofSamples):
    """
    Creating a seperate converting function for Deviations
    as they must be created seperatly
    """
    i = 0 
    j = 0
    shape = np.shape(PlottingDeviation)
    print(shape)
    a = PlottingDeviation.shape[0]
    b = PlottingDeviation.shape[1]
    print("this is from the data conversion function")
    print(a)
    print(b)
    while i < a:
        while j < b:
            begin = NumberofSamples*(i * a + j)
            begin = int(begin) 
            end = NumberofSamples*(i * a + j + 1)
            end = int(end)
            PlottingDeviation[i][j] = hmean(Time[begin: end])
            #the iteration is probably wrong, this is the wrong deviation
            j += 1
        j = 0
        i += 1

    return  PlottingDeviation


#need to fix the samplesize later on
def Prepare(Index, Thread, Time, BestFitnessValue):
    a = len(np.unique(Index))
    b = len(np.unique(Thread))

    #creating several 2 D arrays, which will be needed for the plotting
    PlottingIndex = np.empty((a, b)) # probaly convert this to integer
    PlottingThread =  np.empty((a, b)) # probaly convert this to integer
    PlottingTime =  np.empty((a, b))
    PlottingBestFitnessValue = np.empty((a, b))
    PlottingDeviation = np.empty((a, b))

    print("This is from the preperation function")
    print(PlottingDeviation.shape[0])
    print(PlottingDeviation.shape[1])
    Samplesize = numberofsamplesforIndexandThreads(Index, Thread)    
    #the following equation must hold: NumberofSample * a * b
    #call different converting function to get data in the right format
    ConvertDataDeviations(Time, PlottingDeviation, Samplesize)
    ConvertData(Index, PlottingIndex, Samplesize)
    ConvertData(Thread, PlottingThread, Samplesize)
    ConvertData(Time, PlottingTime, Samplesize)
    ConvertData(BestFitnessValue, PlottingBestFitnessValue, Samplesize)

    Index = PlottingIndex.copy()
    Thread = PlottingThread.copy()
    Time = PlottingTime.copy()
    BestFitnessValue = PlottingBestFitnessValue.copy()
    return PlottingDeviation


def main():
    """
    I radically changed the build up of my function. I propose the new following build up.
    1. A seperate readin function for readability
    2. A Sorting function that call the 2 different sorting patterns at the same time.
    3. A preperation function that resorts the storage of the function so that it will be easier to access for plotting 
    4. A function for plotting with several indicies
    """
    nameoffile = input("What is the location of your file! ")
    #Index, Thread, Time, BestFitnessValue = ReadIn(filelocation)
    Index , Thread , Time , BestFitnessValue = np.loadtxt( nameoffile , delimiter=' ', unpack=True)
    Sorting(Index, Thread, Time, BestFitnessValue)
    Deviations = Prepare(Index, Thread, Time, BestFitnessValue) 



if __name__ == '__main__':
    main()
