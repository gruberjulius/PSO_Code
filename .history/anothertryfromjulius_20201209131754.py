import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def plotJulius(Times, Threads, Deviations, name):
    plt.cla()
    #plt.plot(Times, Threads)
    plt.title(name)
    plt.xlabel("Number of Threads")
    plt.ylabel("Average Times and Standard Deviation needed", rotation = 0)
    plt.savefig(name+'.pdf')




#we need to give the function the number of samples per set
def cleaning(Threads, Times, stepsize):
    size = len(Threads)
    length = size/stepsize 
    #get the number of different threads sizes used
    modifiedTimes = np.empty(0)
    modifiedDeviations = np.empty(0)
    modifiedThreads = np.arange(length) 
    i = 1

    while i < size:
        j = 0
        helparray = np.empty(stepsize)
        if(Threads[i - 1] == Threads[i]):
            helparray[j] = Times[i - 1]
            j += 1
        else:
            np.append(modifiedDeviations , np.std(helparray))
            np.append(modifiedTimes, np.average(helparray))
            j = 0
        i += 1
    
    return modifiedThreads, modifiedTimes, modifiedDeviations

#you still need to add a adding a function name
#you still need to be able the number of samples made
def main():
    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                             "results.dat", delimiter = ' ', unpack = True)


    # function that cleans the data 
    mThreads, mTimes, mDeviations = cleaning(Threads, Time, 10) 

    #plotting function 
    plotJulius(mTimes, mThreads, mDeviations, "SomePlots")


if __name__ == '__main__':
    main()
