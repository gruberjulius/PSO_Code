import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


#Times is the average Time needed
#Threads is the number of threads
#Deviation is the standard deviation
#name is the name of the plot
#xl stands for xlabel to label the x scale
#yl stands for ylabel to label the y scale
def plotJulius(Times, Threads, Deviations, name, xl, yl):
    plt.cla()
    plt.errorbar(Threads,Times,Deviations, linestyle='None', marker='^')
    plt.title(name)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation = 0)
    plt.savefig(name+'.pdf')

#we need to give the function the number of samples per set
#numberoffuncs:number of threads input to the cleaning function i.e.15


#this functions "cleans" the data
#it taked in the big files with number of threads, the times, the standard step size used
#and the number of different functions one wants to plot
#it outputs a cleaned version, hence three arrays of size numoffunc
#the first is the numbering of threads, the second one with the averaged times for each threat step
#and the third one for the standarddeviations for each timestep
def cleaning(Threads, Times, stepsize, numoffunc):
    size = len(Threads)

    #get the number of different threads sizes used
    modifiedTimes = np.empty(numoffunc)
    modifiedDeviations = np.empty(numoffunc)
    modifiedThreads = np.arange(numoffunc) 
    i = 0
    
    while(i < numoffunc):
        modifiedTimes[i] = np.average(Times[(i*stepsize):((i*stepsize)+stepsize)])
        modifiedDeviations[i] = np.std(Times[(i*stepsize):((i*stepsize)+stepsize)])
        i += 1   

    return modifiedThreads, modifiedTimes, modifiedDeviations

#you still need to add a adding a function name
#you still need to be able the number of samples made
#you need to be able to add the name of the data set we use
#stepsize as input variable

def main():
    location = input("The location of the file")

    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                            #"results.dat", delimiter = ' ', unpack = True)
                            location, delimiter = ' ', unpack = True)

    stepsize = input("Number of runs for every thread") #input_a = int(input_a) 
    stepsize = int(stepsize)
    nameoffunc = input("What is the name of the function") 

    #number of runs for every thread
    #stepsize = 10

    #number for how many different types of threads are used
    #should be calculated with numberofthreads = len(Threads)//stepsize 
    #numberofthreads = 15
    numberofthreads = len(Threads)/stepsize
    # function that cleans the data 
    mThreads, mTimes, mDeviations = cleaning(Threads, Time, 
                                    stepsize, numberofthreads) 

    #plotting function 
    plotJulius(mTimes, mThreads, mDeviations,nameoffunc
                ,"Number of Threads","Times in Sec")

if __name__ == '__main__':
    main()
