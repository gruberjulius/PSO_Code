import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
import copy
import sys

#def ReadIn(storagelocation):
    

def plot(Index, Thread, Time, Deviation, storagelocation ='', 
        nameoftheplot='', xl='Chicken Legs', yl = 'Threads'):
    #I will create a seperate plot for every index and put them together in the end
    
    plt.cla()
    plt.figure(figsize=(12,6)) # function to scale the plot how you want it
    #plt.axis([0, 17, 3, 4]) # manually scale the plotting
    plt.title(nameoftheplot)
    plt.grid(linewidth = "1", b=True, which='both', alpha=1)
    plt.xlabel(xl)
    plt.ylabel(yl, rotation=0) 
    plt.savefig( storagelocation + nameoftheplot + '.pdf' )

def main():
    """
    Plotting function
    """
    storagelocation = input("What is the path to your file? " )
    #Index, Thread, Time, Deviation, BestFV = ReadIn(storagelocation)
    Index , Thread , Time , BestFitnessValue = np.loadtxt( nameoffile , delimiter=' ', unpack=True)
    plot(Index, Thread, Time, Deviation) #need to add a the place to store the file

if __name__ == '__main__':
    main()