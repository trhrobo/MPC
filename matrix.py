import numpy as np

a = ([])
print(a)
b = np.array([[0, 1], 
              [2, 3]])
print(b)
a.append(b)
print(a)
c = np.array([[4, 5], 
              [6, 7]])
a.append(c)
print(a)