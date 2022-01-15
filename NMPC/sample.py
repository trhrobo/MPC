import numpy as np
import matplotlib.pyplot as plt
import math

if __name__ == '__main__':
    vs=np.zeros((3, 4))
    vs[0][0]=1
    vs[0][1]=2
    vs[0][2]=3
    vs[0][3]=4
    vs[1][0]=5
    vs[1][1]=6
    vs[1][2]=7
    vs[1][3]=8
    vs[2][0]=9
    vs[2][1]=10
    vs[2][2]=11
    vs[2][3]=12
    print(vs)
    print("----")
    for i in range(3):
        print(i)
        print("num")
        print(vs[:, :i-1])    
        print("----")