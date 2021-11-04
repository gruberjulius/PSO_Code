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
import plottingwithdeviationfile
import plottingclassicalyfile
import plottingwithfitnessfile
import buildupfunctions


#PlottingIndex1, PlottingThread1, PlottingTime1, PlottingBestFitnessValue1, PlottingDeviation1 =\
#                    Deviations = Prepare(Index1, Thread1, Time1, BestFitnessValue1) 


class dataclass():

    def __init__(self, storagelocation):
        self._Index, self._Thread, self._Time, self._BestFV = np.loadtxt( storagelocation , delimiter=' ' , unpack=True)

    def readin(storagelocation,self):
        self._Index, self._Thread, self._Time, self._BestFV = np.loadtxt( storagelocation , delimiter=' ' , unpack=True)


    def ConvertDataDeviations(self, NumberofSamples):
    i = 0 
    j = 0
    shape = np.shape(self._PlottingDeviation)
    a = self._PlottingDeviation.shape[0]
    b = self._PlottingDeviation.shape[1]
    while i < a:
        while j < b:
            begin = NumberofSamples*(i * b + j)
            begin = int(begin) 
            end = NumberofSamples*(i * b + j + 1)
            end = int(end)
            self._PlottingDeviation[i][j] = np.std(self._Time[begin: end])
            #the iteration is probably wrong, this is the wrong deviation
            j += 1
        j = 0
        i += 1

    #def prepare(self._Index, self._Thread, self._Time, self._BestFV):
    def prepare(self):
        a = len(np.unique(self._Index))
        b = len(np.unique(self._Thread))

        #creating several 2 D arrays, which will be needed for the plotting
        PlottingIndex = np.empty((a, b)) # probaly convert this to integer
        PlottingThread =  np.empty((a, b)) # probaly convert this to integer
        PlottingTime =  np.empty((a, b))
        PlottingBestFitnessValue = np.empty((a, b))
        PlottingDeviation = np.empty((a, b))
        Samplesize = 50 
        #####number of samples and threads should be here    
        #the following equation must hold: NumberofSample * a * b
        #call different converting function to get data in the right format
        ConvertDataDeviations(self, Samplesize)
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

    if typeofplot== c:
        #call the original plotting function
        #you already have this one
        plottingwithdeviation(Datasets)

    if typeofplot== d:
        #call the deviation function
        plottingwithdeviation(Datasets)

    if typeofplot== f:
        #calls the fitness plotting function
        plottingwithbestfitness(Datasets)
            

if __name__ == '__main__':
    main()
