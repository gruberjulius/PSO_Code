import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
import copy
import sys
import buildupfunctions


def plottingwithdeviation(Datasets):
    numberofdatasets = len(Datasets)
    typeofplot = input("What kind of plot do you want to do, c for classical, d for deviation, f for fitness")
    nameofplot = input("What is the name of the plot! ")
    storagelocation = input("What is the location of the storage? ")
    xlabel = input("What is the label for the x axis? ")
    ylabel = input("What is the label for the y axis? ")
    titleofplot = input("What is the title of the plot? ") 
    startindex = int(input("What is the starting index of the functions you want to plot: "))
    endindex = int(input("what is the end index of functions U use? ")) # exklusiv

    #copy paste the plotting function for the fitness value

    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 35, 0.01, 0.05]) # manually scale the plotting
    plt.title(titleoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=90)

    while j:= 0 < numberofdatasets:
            i = int(startindex)
            k = int(endindex)
            labelname = input("How do you want to label the dataset? ")
            while i < k :#numberofdifferentindicies
                ThreadDummy = Datasets[j]._PlottingThread[i:(i+ 1)]
                TimeDummy = Datasets[j]._PlottingTime[i:(i + 1)]
                Time2Dummy = TimeDummy[0]
                Thread2Dummy = ThreadDummy[0]
                plt.plot(Thread2Dummy, Time2Dummy, label = labelname)
                i += 1
        j+= 1

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + '.pdf' )
