import numpy as np

def countDistinct(arr):

    """
    Takes an numpy array as an input 
    Ouput will be the number of different elements and the different elements in a numpy array
    Tested it works fine
    """  
    res = 0
    n = len(arr)
    # Pick all elements one by one 
    for i in range(1, n): 
        j = 0
        for j in range(i): 
            if (arr[i] == arr[j]): 
                break
  
        # If not printed earlier, then print it 
        if (i == j + 1): 
            res += 1
    return res


def main():
    """
    I radically changed the build up of my function. I propose the new following build up.
    1. A seperate readin function for readability
    2. A Sorting function that call the 2 different sorting patterns at the same time.
    3. A preperation function that resorts the storage of the function so that it will be easier to access for plotting 
    4. A function for plotting with several indicies
    """
    filelocation = input("What is the location of your file! ")
    #Index, Thread, Time, BestFitnessValue = ReadIn(filelocation)
    Index , Thread , Time , BestFitnessValue = np.loadtxt( filelocation , delimiter=' ' , unpack=True)
    Sorting(Index, Thread, Time, BestFitnessValue)
    print(Index)
    print(Thread)


if __name__ == '__main__':
    main()
