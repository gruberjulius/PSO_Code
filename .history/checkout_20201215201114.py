import numpy as np


a = np.array([1,2,3])
b = np.array([4,5,6])

print(a)
print(b)
a[1] = b[2]
print(a)
print(b)
b[0] = a[2]
print(a)
print(b)