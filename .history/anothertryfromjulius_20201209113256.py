import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

#we need to give the function the number of samples per set
def plotting(Threads, Times, name, numberofsamplesperstep):
    length = len(Threads)/numberofsamplesperstep 
    #get the number of iterations we need
    modifiedTimes = np.empty(0)
    modifiedDeviations = np.empty(0)
    i = 1

    for i in len(Threads):
        j = 0
        helparray = np.empty(numberofsamplesperstep)
        if(Threads[i - 1] == Threads[i]):
            np.append(helparray[j], Times[i - 1])
            j +=1
        else:
            np.append(modifiedDeviation , np.std(helparray))
            np.append(modifiedTimes, np.average(helparray))
            j = 0
        i +=1
    
    plt.cla()
    plt.title(name)
    plt.savefig(name+'.pdf')

def main():
    name = "SomePlots"
    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                             "results.dat", delimiter = ' ', unpack = True)
    plotting(Threads, Time, name)
    size = 10

if __name__ == '__main__':
    main()
