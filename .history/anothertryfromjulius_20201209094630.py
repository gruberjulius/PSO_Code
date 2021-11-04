import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

sns.set_theme() #set standard theme for seaborn

def confidence_plot(Threads ,Time, name):
    plt.cla()
    sns.regplot(Threads, Time)
    plt.title(name)
    plt.savefig(name+'.pdf')
    plt.show()

def main():
    name = "SomePlots"
    Index, Threads, Time, BestFitnessValue = np.loadtxt(
                            # "results.dat", delimiter = ' ', unpack = True)
    confidence_plot(Threads, Time, name)


#if __name__ == '__main__':
    main()#