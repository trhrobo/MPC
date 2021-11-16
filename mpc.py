#-*- coding: utf-8 -*-
import numpy as np
import matplotlib.pyplot as plt
import cvxpy

print("Simulation Start")

goal_x = np.array([[0], [0]])
start_x = np.array([[10], [10]])
x = start_x

A = np.array([[1, 1], [0, 1]])
B = np.array([[0.5], [1]])

Q = np.array([[1, 1], [1, 1]])
R = np.array([[1, 1], [1, 1]])

N = 10

H = []
for i in range(N):
    temp_A = np.linalg.matrix_power(A, i)
    H.append(temp_A)
K = []
#行
for i in range(N):
    #列
    val = []
    for t in range(N):
        temp_val_A = np.linalg.matrix_power(A, i)
        if(i >= t):
            temp_val = temp_val_A * B
        else:
            temp_val = temp_val_A
        val.append(temp_val)

    K.append(val)

I = np.eye(2)
bmQ = np.array([[Q, I, I], [I, Q, I], [I, I, Q]])
bmR = np.array([[R, I, I], [I, R, I], [I, I, R]])

print(bmQ)

Phi = 2 * (np.transpose(K) * bmQ * K + bmR)
print(Phi)
Gamma = 2 * np.transpose(x) * np.transpose(H) * bmQ * K

U = cvxpy.Variable(4)
objective = cvxpy.sum_squares(0.5 * U.transpose * Phi * U + Gamma * U)
prob = cvxpy.Problem(cvxpy.Minimize(objective))
prob.solve()
