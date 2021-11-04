import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches

#calculte the harmonic mean of an array
def hmean(arr):
    """calculates the harmonic mean of a functions
    the elements of the array have to be positive in order to get a accurate results
    """
    inverse = 1./arr #invert the array elements
    lenarr = len(inverse) # get the length of the array
    sumarr = sum(inverse) # get the sum of all the inverted elements
    if sumarr != 0 :
        return lenarr/sumarr # return the harmonic mean

#calculate the number of distinct elements in an array
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

#function to calculate the step size for 2 threads and arrays
def numberofsamplesforIndexandThreads(Index, Thread):
    """
    function that returns the number of samples for each index and thread
    Input : Indexlist and Threadlist
    Output: Number of samples per Index and Threas
    """
    checkvariable = 0
    firstIndex = Index[0]
    firstThread = Thread[0]
    i = 0
    while i < len(Index):
        if firstIndex ==  Index[i] and firstThread == Thread[i]:
            checkvariable += 1
        else :
            break
        i += 1

    return checkvariable

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

#function that readsin elements without all the further mess
def ReadIn(filelocation):
    """ 
    This is a readin function for our files. 
    I seperated it from the main file in order to make it more readable
    """
    Index, Thread, Time, BestFitnessValue = np.loadtxt(
        filelocation, delimiter=' ', unpack=True)
    return Index, Thread, Time, BestFitnessValue

#function that calls two sorting functions seperatly
def Sorting(Index, Thread, Time, BestFitnessValue):
    """
    function that calls the 2 sorting functions seperatly. Excluded this for readabilits
    """
    modifiedsortIndex(Index, Thread, Time, BestFitnessValue)
    modifiedsortIndexThreads(Index, Thread, Time, BestFitnessValue)

#function that converts data in out favourite format
#the iteration in there is probably completly wrong
def ConvertData(Format, PlottingFormat, NumberofSamples):
    """
    Converts Data to different format to make it easier to plot
    """
    i = 0
    j = 0
    a = np.size(PlottingFormat, axis = 0)
    b = np.size(PlottingFormat, axis = 1)
    while i < a:
        while j < b:
            PlottingFormat[i][j] = hmean(Format[(NumberofSamples*(i + j)):(NumberofSamples*(i + j + 1))])
            #the iteration is probably wrong
            j += 1
        i += 1
    
    Format = PlottingFormat.copy()


def ConvertDataDeviations(Time, PlottingDeviation, NumberofSamples):
    """
    Creating a seperate converting function for Deviations
    """
    i = 0 
    j = 0
    a = np.size(PlottingDeviation, axis = 0)
    b = np.size(PlottingDeviation, axis = 1)
    while i < a:
        while j < b:
            PlottingDeviation[i][j] = np.std(Time[(i + j):(i + j + 1)])
            #the iteration is probably wrong
            j += 1
        i += 1

    return  PlottingDeviation


def Prepare(Index, Thread, Time, BestFitnessValue):
    a = countDistinct(Index)
    b = countDistinct(Thread)


    #creating several 2 D arrays, which will be needed for the plotting

    PlottingIndex = np.array((a, b)) # probaly convert this to integer
    PlottingThread =  np.array((a, b)) # probaly convert this to integer
    PlottingTime =  np.array((a, b))
    PlottingBestFitnessValue = np.array((a, b))
    PlottingDeviation = np.array(((a, b))

    NumberofSamples = numberofsamplesforIndexandThreads(Index, Thread) 
    
    #the following equation must hold: NumberofSample * a * b
    #call different converting function to get data in the right format
    ConvertDataDeviations(Time, PlottingDeviation, NumberofSamples)
    ConvertData(Index, PlottingIndex, NumberofSamples)
    ConvertData(Thread, PlottingThread, NumberofSamples)
    ConvertData(Time, PlottingTime, NumberofSamples)
    ConvertData(BestFitnessValue, PlottingBestFitnessValue, NumberofSamples)

    Index = PlottingIndex.copy()
    Thread = PlottingThread.copy()
    Time = PlottingTime.copy()
    BestFitnessValue = PlottingBestFitnessValue.copy()
    return PlottingDeviation


"""def Plot(Index, Thread, Time, BestFitnessValue, Deviations, storagelocation):
    # you have to create subplots for every element
    Deviations_Plot = mpatches.Patch(color ='gray', label ='Deviations')
    Median_Plot = mpatches.Patch(color ='black', label='Median')
    plt.cla()

    for i in range(Index)
"""



def main():
    """
    I radically changed the build up of my function. I propose the new following build up.
    1. A seperate readin function for readability
    2. A Sorting function that call the 2 different sorting patterns at the same time.
    3. A preperation function that resorts the storage of the function so that it will be easier to access for plotting 
    4. A function for plotting with several indicies
    """
    filelocation = input("What is the location of your file! ")
    Index, Thread, Time, BestFitnessValue = ReadIn(filelocation)
    Sorting(Index, Thread, Time, BestFitnessValue)
    #function that prepares the elements for plotting & creates a deviation file
    Deviations = Prepare(Index, Thread, Time, BestFitnessValue) 
    storagelocation = input("Input the storage location, where you want to store your files? ")
    Plot(Index, Thread, Time, BestFitnessValue, Deviations, storagelocation)




if __name__ == '__main__':
    main()
