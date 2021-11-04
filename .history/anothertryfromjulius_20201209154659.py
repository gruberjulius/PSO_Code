import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#xl stands for xlabel to label the x scale
#yl stands for ylabel to label the y scale
def plotJulius(Times, Threads, Deviations, name, xl, yl):
    plt.cla()
    #plt.plot(Times, Threads)
    plt.errorbar(Threads,Times,Deviations, linestyle='None', marker='^')
    plt.title(name)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation = 0)
    plt.savefig(name+'.pdf')




#we need to give the function the number of samples per set
#numberoffuncs:number of threads input to the cleaning function i.e.15
def cleaning(Threads, Times, stepsize, numoffunc):
    size = len(Threads)
    
    #get the number of different threads sizes used
    modifiedTimes = np.empty(numoffunc)
    modifiedDeviations = np.empty(numoffunc)
    modifiedThreads = np.arange(numoffunc) 
    i = 0

    
    """
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
    """
    
    while(i < numoffunc):
        modifiedTimes[i] = np.average(Times[(i*10):((i*10)+10)])
        modifiedDeviations[i] = np.std(Times[(i*10):((i*10)+10)])
        i += 1   

    

    return modifiedThreads, modifiedTimes, modifiedDeviations

#you still need to add a adding a function name
#you still need to be able the number of samples made
def main():
    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                             "results.dat", delimiter = ' ', unpack = True)

    #number of samples for every thread
    numberofsamplesforiteration = 10

    #number for how many different types of threads are used
    numberofthreads = 15
    # function that cleans the data 
    mThreads, mTimes, mDeviations = cleaning(Threads, Time, 
                                    numberofsamplesforiteration, numberofthreads) 

    #plotting function 
    plotJulius(mTimes, mThreads, mDeviations, "Plotting Threads against Time needed"
                ,"Number of Threads","y-axes")


if __name__ == '__main__':
    main()
