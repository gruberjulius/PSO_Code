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


def plottingwithclassical(Datasets):
    numberofdatasets = len(Datasets)
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
    plt.title(titleofplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel, rotation=90)
    j = 0
    while j < numberofdatasets:
            i = int(startindex)
            k = int(endindex)
            l = str(j)
            labelname = input("How do you want to label the dataset" + l + "? ")
            while i < k :#numberofdifferentindicies
                ThreadDummy = Datasets[j]._PlottingThread[i:(i+ 1)]
                print("Beginning of printing of data set")
                print(ThreadDummy)
                TimeDummy = Datasets[j]._PlottingTime[i:(i + 1)]
                print(TimeDummy)
                Time2Dummy = TimeDummy[0]
                print(Time2Dummy)
                Thread2Dummy = ThreadDummy[0]
                print(Thread2Dummy)
                plt.plot(Thread2Dummy, Time2Dummy, label = labelname)
                i += 1
            j+= 1

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + '.pdf' )
