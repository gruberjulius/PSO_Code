# general idea
# main function that calls the needed classes
# data class, where the dat is stpred
#  
# data class has a read in function to read in the data
# data class it stores the data as a member, seperate for every data set
# plotting class is a befreiended class of the function
# plotting class has a function that plots the functions in a normal plot
# plotting class has a function that plots the function with a deviaton


#you need to copy paste the prepare function

#you need to check concerning string modification it would make things much easier

#put this into a seperate file, which you can call then
# put each plotting function into a seperate file

import plottingwithdeviationfile
import plottingclassicaly
import plottingwithfitnessfile

class dataclass:
    #you need to add the conversion to deviations automatically

    def __init__(self):
        self



    def readin(Dataset):
        self._Index, self._Thread, self.Time, self.BestFV = 
            np.loadtxt(Dataset,delimiter=' ', unpack=True)
    

class plottingwithdeviation:



class plottingwithbestfitness:



class plottingclassicclass:

    def init(self):

    def plot1set(Dataset):
        #fill in the material for 1 dataset copy paste the old one

def constfunction(constant value):
    """
    function that gets as input  a constant value
    and plots the value on the whole axes
    """


    

def main():
    numberofdatasets = int(input("How many functions do you want to plot! "))

    #create a while loop that files n data sets and calls them accoordinggly
    #create a seperate counting variable later on, 
    #you will need the one again to iterate again over the functions
    #still need to manually add the file paths
    
    #something like this should work out
    while i:=0 < numberofdatasets:
        dataclass_ += "i"
        dataclass(dataclass_) # convert to type data class
        locationofdataset = input("What is the location of the data set? ")
        dataclass_ + "i" = dataclass.readin(locationofdataset)
        i+=1
    
    typeofplot = input("What kind of plot do you want to do, c for classical, d for deviation, f for fitness")
    nameofplot = input("What is the name of the plot! ")
    storagelocation = input("What is the location of the storage? ")
    xlabel = input("What is the label for the x axis? ")
    ylabel = input("What is the label for the y axis? ")
    titleofplot = input("What is the title of the plot? ") 
    startindex = int(input("What is the starting index of the functions you want to plot: "))
    while j:= 0 < numberofdatasets:
        if typeofplot==c:
            #call the original plotting function
            #you already have this one
            plottingwithdeviation(numberofdatasets)

        if typeofplot==d:
            #call the deviation function
            plottingwithdeviation(numberofdatasets)

        if typeofplot==f:
            #calls the fitness plotting function
            plottingwithbestfitness(numberofdatasets)
            

if __name__ == '__main__':
    main()
