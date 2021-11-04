# general idea
# main function that calls the needed classes
# data class, where the dat is stpred
#  
# data class has a read in function to read in the data
# data class it stores the data as a member, seperate for every data set
# plotting class is a befreiended class of the function
# plotting class has a function that plots the functions in a normal plot
# plotting class has a function that plots the function with a deviaton


#you need to copy paste the prepare function

#you need to check concerning string modification it would make things much easier

#put this into a seperate file, which you can call then
# put each plotting function into a seperate file
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
import copy
import sys
from plottingwithdeviationfile import *#plottingwithdeviation
from plottingclassicalyfile import *
from plottingwithfitnessfile import *
from  buildupfunctions import *


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
            PlottingFormat[i][j] =  hmean(Format[50 * (i  * b + j) : (i * b + j + 1)*50])
            j += 1
        j = 0
        i += 1

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

class dataclass():

    def __init__(self, storagelocation):
        self._Index, self._Thread, self._Time, self._BestFV = np.loadtxt( storagelocation , delimiter=' ' , unpack=True)


    #def prepare(self._Index, self._Thread, self._Time, self._BestFV):
    def prepare(self):
        a = len(np.unique(self._Index))
        b = len(np.unique(self._Thread))

        #creating several 2 D arrays, which will be needed for the plotting
        self._PlottingIndex = np.empty((a, b)) # probaly convert this to integer
        self._PlottingThread =  np.empty((a, b)) # probaly convert this to integer
        self._PlottingTime =  np.empty((a, b))
        self._PlottingBestFitnessValue = np.empty((a, b))
        self._PlottingDeviation = np.empty((a, b))
        Samplesize = 50 
        #####number of samples and threads should be here    
        #the following equation must hold: NumberofSample * a * b
        #call different converting function to get data in the right format
        ConvertDataDeviations(self._Time, self._PlottingDeviation, Samplesize)
        ConvertData(self._Index, self._PlottingIndex)
        ConvertData(self._Thread, self._PlottingThread)
        ConvertData(self._Time, self._PlottingTime)
        ConvertData(self._BestFV, self._PlottingBestFitnessValue)
        #we return all elements of the computation seperatly in order to call out plotting function later
        return self._PlottingIndex, self._PlottingThread, self._PlottingTime, self._PlottingBestFitnessValue, self._PlottingDeviation



    


def constfunction(constant_value):
    """
    function that gets as input  a constant value
    and plots the value on the whole axes
    """


    

def main():
    numberofdatasets = int(input("How many functions do you want to plot! "))

    #create a while loop that files n data sets and calls them accoordinggly
    #create a seperate counting variable later on, 
    #you will need the one again to iterate again over the functions
    #still need to manually add the file paths
    
    #something like this should work out

    #create a wonderful list of the data sets
    Datasets = []
    
    i = 0
    while i  < numberofdatasets:
        locationofdataset = input("What is the location of the data set? ")
        dummy = dataclass(locationofdataset)
        #dummy.readin(locationofdataset)
        Datasets.append(dummy)
        Datasets[len(Datasets) - 1].prepare()
        i += 1
    
    typeofplot = input("What kind of plot do you want to do: c for classical, d for deviation, f for fitness")

    if typeofplot== "c":
        #call the original plotting function
        #you already have this one
        plottingwithclassical(Datasets)

    if typeofplot== "d":
        #call the deviation function
        plottingwithdeviation(Datasets)

    if typeofplot== "f":
        #calls the fitness plotting function
        plottingwithbestfitness(Datasets)
            

if __name__ == '__main__':
    main()
