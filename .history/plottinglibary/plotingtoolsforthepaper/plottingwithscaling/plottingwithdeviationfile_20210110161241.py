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
    print("a")
    nameofplot = input("What is the name of the plot! ")
    storagelocation = input("What is the location of the storage? ")
    xlabel = input("What is the label for the x axis? ")
    ylabel = input("What is the label for the y axis? ")
    titleofplot = input("What is the title of the plot? ") 
    startindex = int(input("What is the starting index of the functions you want to plot: "))
    endindex = int(input("what is the end index of functions U use? ")) # exklusiv
    #print("b")
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    lower_bound = input("What is the lower bound?")
    lower_bound = float(lower_bound)
    upper_bound = input("What is the upper bound?")
    upper_bound = float(upper_bound)
    lower_height = input("What is the lower height?")
    lower_height = float(lower_height)
    upper_height = input("What is the upper heigh?")
    upper_height = float(upper_height)
    plt.axis([lower_bound, upper_bound, lower_height, upper_height])
    #plt.axis([0, 34, 0, 1.75]) # manually scale the plotting
    plt.title(titleofplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel, rotation=90) 
    #print("c")

    j = 0
    while j < numberofdatasets:
            i = int(startindex)
            k = int(endindex)
            l = str(j)
            labelname = input("How do you want to label the dataset" + l + "? ")
            constant = input("is it contant? y for yes, n for no! ")
            if constant =="y":
                ThreadDummy = Datasets[j]._PlottingThread[i:(i+ 1)]
                #print("Beginning of printing of data set")
                #print(ThreadDummy)
                TimeDummy = Datasets[j]._PlottingTime[i:(i + 1)]
                #print(TimeDummy)
                Time2Dummy = TimeDummy[0]
                #print(Time2Dummy)
                Thread2Dummy = ThreadDummy[0]
                tmp_time = Time2Dummy[0]
                tmp_thread = int(Thread2Dummy[0])
                tmp_thread_size = 33
                nTime2Dummy = np.full((tmp_thread_size), tmp_time)
                nThread2Dummy = np.zeros(tmp_thread_size)
                for z in range(tmp_thread_size):
                    nThread2Dummy[z] = z + 1
                plt.plot(nThread2Dummy, nTime2Dummy, label = labelname)#, fmt='o')
                i += 1
            elif constant == "n":
                while i < k :#numberofdifferentindicies
                    ThreadDummy = Datasets[j]._PlottingThread[i:(i+ 1)]
                    TimeDummy = Datasets[j]._PlottingTime[i:(i + 1)]
                    DeviationDummy = Datasets[j]._PlottingDeviation[i:(i+ 1)]
                    #print(DeviationDummy)
                    Deviation2Dummy = DeviationDummy[0]
                    #print("AAAAAAAAAAAAAAAAAAA")
                    #print(Deviation2Dummy)
                    Time2Dummy = TimeDummy[0]
                    Thread2Dummy = ThreadDummy[0]
                    plt.errorbar(Thread2Dummy, Time2Dummy, Deviation2Dummy, label = labelname, fmt='o', elinewidth=1.5, capsize=5, capthick=3)
                    #, uplims=upperlimits, lolims = lowerlimits)#, barsabove = 'True')#, uplims=True, lolims=True)#, linestyle='None', fmt='o',  label = labelname)
                    i += 1
            j += 1
    
    #print("d")

    plt.legend(loc='upper right')#, handles = [Median_Plot, Deviations_Plot])
    plt.savefig( storagelocation + '.pdf' )