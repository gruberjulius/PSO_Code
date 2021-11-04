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
import plottingclassicaly
import plottingwithfitnessfile



class dataclass(obejcts):
    #you need to add the conversion to deviations automatically


    def readin(Dataset):
        self._Index, self._Thread, self.Time, self.BestFV = np.loadtxt(Dataset,delimiter=' ', unpack=True)
    


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
    

    while i := 0 < numberofdatasets:
        locationofdataset = input("What is the location of the data set? ")
        dataclass dummy  
        dummy.readin(locationofdataset)
        Datasets.append(dummy)
        i += 1
    
    typeofplot = input("What kind of plot do you want to do: c for classical, d for deviation, f for fitness")

    if typeofplot==c:
        #call the original plotting function
        #you already have this one
        plottingwithdeviation(Datasets)

    if typeofplot==d:
        #call the deviation function
        plottingwithdeviation(Datasets)

    if typeofplot==f:
        #calls the fitness plotting function
        plottingwithbestfitness(Datasets)
            

if __name__ == '__main__':
    main()
