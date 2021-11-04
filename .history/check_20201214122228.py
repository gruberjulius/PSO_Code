import numpy as np

arr = np.array(

def hmean(arr):
    
    inverse = 1./arr #invert the array elements
    lenarr = len(inverse) # get the length of the array
    sumarr = sum(arr) # get the sum of all the inverted elements
    if sumarr != 0 :
        return lenarr/sumarr # return the harmonic mean
