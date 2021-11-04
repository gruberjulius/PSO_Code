import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#we need to give the function the number of samples per set
def plotting(Threads, Times, name, numberofsamplesperstep):
    length = len(Threads)/numberofsamplesperstep 
    #get the number of different threads sizes used
    modifiedTimes = np.empty(0)
    modifiedDeviations = np.empty(0)
    modifiedThreads = np.arange(length) 
    i = 1

    while i < len(Threads):
        j = 0
        helparray = np.empty(numberofsamplesperstep)
        if(Threads[i - 1] == Threads[i]):
            np.append(helparray[j], Times[i - 1])
            j += 1
        else:
            np.append(modifiedDeviations , np.std(helparray))
            np.append(modifiedTimes, np.average(helparray))
            j = 0
        i += 1
        
    plt.cla()
    plt.plot(modifiedTimes, modifiedThreads)
    plt.title(name)

    plt.xlabel("Number of Threads")
    plt.ylabel("Average Times and Standard Deviation needed", rotation = 0)
    plt.savefig(name+'.pdf')

#you still need to add a adding function 
#you still need to be able the number of samples made
def main():
    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                             "results.dat", delimiter = ' ', unpack = True)
    plotting(Threads, Time, "SomePlots", 10)

if __name__ == '__main__':
    main()
