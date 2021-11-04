import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
#import pandas as pd
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import matplotlib.patches as mpatches
import copy
import sys

def plottingwithbestfitness(Datasets):
    numberofdatasets = len(Datasets)
    typeofplot = input("What kind of plot do you want to do, c for classical, d for deviation, f for fitness")
    nameofplot = input("What is the name of the plot! ")
    storagelocation = input("What is the location of the storage? ")
    xlabel = input("What is the label for the x axis? ")
    ylabel = input("What is the label for the y axis? ")
    titleofplot = input("What is the title of the plot? ") 
    startindex = int(input("What is the starting index of the functions you want to plot: "))
    endindex = int(input("what is the end index of functions U use? ")) # exklusiv

    while j:= 0 < numberofdatasets:
        #copy paste the plotting routing for the fitness value


        j+= 1