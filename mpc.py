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

#H設定
#サイズの確認
H = A
Adi = A
for i in range(N):
    Adi = np.dot(Adi, A)
    H = np.append(H, Adi)
#サイズの確認
H = np.reshape(H, (2*N, 2))
    #temp_A = np.linalg.matrix_power(A, i)
    #H.append(temp_A)

K = []
for i in range(N):
    temp_val = [B, B, B, B, B, B, B, B, B, B]
    for t in range(N):
        temp_val_A = np.linalg.matrix_power(A, i)
        if(i == t):
            temp_val[t] = B
        elif(i > t):
            temp_val[t] = np.dot(temp_val_A, B)
        else:
            zero = np.array([[0], [0]])
            temp_val[t] = zero
   
    val = np.array([[temp_val[0], temp_val[1], temp_val[2], temp_val[3], temp_val[4], temp_val[5], temp_val[6], temp_val[7], temp_val[8], temp_val[9]]])
    #print(val)
    K.append(val)

#サイズ確認する
#二次計画問題
zero = np.zeros(2)
bmQ = np.zeros((2*N, 2*N))
bmR = np.zeros((N, N))
for i in range(2*N):
    bmQ[i][i] = 1
for i in range(N):
    bmR[i][i] = 1
    
Phi = 2 * (np.dot(np.transpose(K) * np.dot(bmQ * K)) + bmR)
Gamma = 2 * np.dot(np.transpose(x) * np.dot(np.transpose(H) * np.dot(bmQ * K)))


Time = 100
for i in range(Time):
    #拘束条件無し
    #U = cvxpy.Variable(4)
    objective = cvxpy.Minimize(0.5*cvxpy.quad_form(x, P)+gm.T@x)
    problem = cvxpy.Problem(objective)
    sol = problem.solve()    
    
#グラフ表示

