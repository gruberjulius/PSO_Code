import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def plotting(Threads, Times,name):
    modifiedTimes
    modifiedDeviations
    i = 1
    for i in len(Threads):
        helparray
        if(Threads[i - 1] == Threads[i]):
            np.append(helparray, Threads[i - 1])
        else:
            np.append(modifiedDeviations, np.std(helparray))
            np.append(modifiedTimes, np.average(helparray))
            helparray
        i +=1
    
    plt.cla()

    plt.title(name)
    plt.savefig(name+'.pdf')
    


def confidence_plot(Threads ,Time, name):
    plt.cla()
    
    plt.title(name)
    plt.savefig(name+'.pdf')
    #plt.show()

def main():
    name = "SomePlots"
    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                             "results.dat", delimiter = ' ', unpack = True)
    confidence_plot(Threads, Time, name)

if __name__ == '__main__':
    main()
