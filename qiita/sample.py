import numpy as np
import matplotlib.pyplot as plt
import math
import time

if __name__ == '__main__':
    vs=np.zeros((3,3))
    vs[0][0]=1;
    vs[0][1]=2;
    vs[0][2]=3;
    vs[1][0]=4;
    vs[1][1]=5;
    vs[1][2]=6;
    vs[2][0]=7;
    vs[2][1]=8;
    vs[2][2]=9;
    print(vs)
    print(vs[:, :2])