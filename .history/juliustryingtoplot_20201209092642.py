import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
#import csv
import pandas as pd
import seaborn as sns
#sns.set_theme(style="darkgrid")

def plot_2D(N, R, name,xl,yl):
    plt.cla()
    plt.plot(N, R)
    #plt.plot(N, np.sqrt(N)) #auskommentiert für S
    plt.xlabel(xl)
    plt.ylabel(yl, rotation = 0)
    #plt.legend(["MD","sqrt(N)"]) #auskom für S
    plt.title(name)
    plt.savefig(name + '.pdf')

def calculateaverage(Time):
    return sum(Time)/len(Time)

def confidence_plot_thread_vs_time(Threads, Time, name):
    sns.regplot(Threads, Time)
    sns.relplot(x="Threads", y= "Time needed", 

#def main():
#    Index, T, R, S = np.loadtxt('results.dat', delimiter = ' ', unpack = True)
#    return 1

def main():
    print("Hello World!")

#def  main(hello):
    #you have to add input output so that we can read in the functions name
 #   name = "Sample Function"
    #r0 = 0.00154
    #N, T, R, S = np.loadtxt('results.dat', delimiter = ' ', unpack = True)
    #
   # plot_2D(S, R, "end2end_R2S","Steps","R(S)")
  #  Index, Threads, Time, BestFitnessValue = np.loadtxt(
                            # "results.dat", delimiter = ' ', unpack = True)
    #Time = np.loadtxt('results.dat', delimiter = ' ')
    #plot_2D(Threads, Time, 
    #"Plotting number of Thread against Time", "Number of Threads", "Time")

   # confidence_plot_thread_vs_time(Threads, Time, name) 

#def main():
 #   name = "Sample Function"
#
#    Index, Threads, Time, BestFitnessValue = np.loadtxt("results.dat", delimiter = ' ', unpack = True)

 #   confidence_plot_thread_vs_time(Threads, Time, name)


if __name__ == '__main__':
    main()
