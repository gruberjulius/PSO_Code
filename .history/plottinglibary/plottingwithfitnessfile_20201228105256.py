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

def plottingwithbestfitness(Datasets):
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

            labelname = input("How do you want to label the dataset? ")
            
            constant = input("is it contant? y for yes, n for no! ")
            if constant =="y":
                ThreadDummy = Datasets[j]._PlottingThread[i:(i+ 1)]
                #print("Beginning of printing of data set")
                #print(ThreadDummy)
                BestFVDummy = Datasets[j]._PlottingBestFitnessValue[i:(i + 1)]
                #print(TimeDummy)
                BestFV2Dummy = BestFVDummy[0]
                #print(Time2Dummy)
                Thread2Dummy = ThreadDummy[0]
                tmp_FV = BestFV2Dummy[0]
                tmp_thread = int(Thread2Dummy[0])
                tmp_thread_size = 33
                nBestFV2Dummy = np.full((tmp_thread_size), tmp_FV)
                nThread2Dummy = np.zeros(tmp_thread_size)
                for z in range(tmp_thread_size):
                    nThread2Dummy[z] = z + 1
                plt.plot(nThread2Dummy, nBestFV2Dummy, label = labelname)
                i += 1
            elif constant == "n":
            
                while i < k :#numberofdifferentindicies
                    ThreadDummy = Datasets[j]._PlottingThread[i:(i+ 1)]
                    BestFVDummy = Datasets[j]._PlottingBestFitnessValue[i:(i + 1)]
                    BestFV2Dummy = BestFVDummy[0]
                    Thread2Dummy = ThreadDummy[0]
                    plt.plot(Thread2Dummy, BestFV2Dummy, label = labelname)
                    i += 1
            j+= 1

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + '.pdf' )