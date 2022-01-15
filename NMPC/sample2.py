import numpy as np
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    e=np.zeros((3,4))
    a=np.zeros((3,1))
    e[0, 0]=1
    e[0, 1]=2
    e[0, 2]=3
    e[0, 3]=4
    e[1, 0]=5
    e[1, 1]=6
    e[1, 2]=7
    e[1, 3]=8
    e[2, 0]=9
    e[2, 1]=10
    e[2, 2]=11
    e[2, 3]=12
    a[0, 0]=1
    a[1, 0]=2
    a[2, 0]=3
    print(e)
    print(np.linalg.pinv(e))
    print(np.dot(e,np.linalg.pinv(e)))