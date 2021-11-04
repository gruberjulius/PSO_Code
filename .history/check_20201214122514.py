import numpy as np

arr = np.array([1,4,4])


def hmean(arr):
    print(arr)
    inverse = 1./arr #invert the array elements
    print(arr)
    lenarr = len(inverse) # get the length of the array
    print(lenarr)
    sumarr = sum(arr) # get the sum of all the inverted elements
    if sumarr != 0 :
        return lenarr/sumarr # return the harmonic mean


bro = hmean(arr)

print(bro)