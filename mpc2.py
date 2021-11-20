#-*- coding: utf-8 -*-
import numpy as np
import scipy
from scipy import signal 
import cvxpy
import matplotlib.pyplot as plt
c = 1
k = 1
A = np.array([[-c, k], [1, 0]])
B = np.array([[1], [0]])
C = np.array([0, 1])
D = np.array([1])

#離散化
dt = 0.01
dsys = scipy.signal.cont2discrete((A, B, C, D), dt, method = 'zoh')
Ad = dsys[0]
Bd = dsys[1]
Cd = dsys[2]

print(Ad)
print(Bd)
print(Cd)

N = 10

#各種定数の設定
#Phi
Adi = Ad
Phi = Ad
for i in range(N-1):
    Adi = np.dot(Adi, Ad)
    Phi = np.append(Phi, Adi)
    
Phi = np.reshape(Phi, (2*N,2))

#Ganma
Ganma = Bd
Adi = Ad
for i in range(N - 1):
    Ganma = np.append(Ganma, np.dot(Adi, Bd))
    Adi = np.dot(Adi, Ad)
    
Ganma = np.reshape(Ganma, (2*N,1))

#Theta
#分かってない
zero = np.zeros((2,1))
Theta = Ganma
Ganma2 = Ganma 
for i in range(N-1):
    Ganma2 = np.append(zero, Ganma2, 0)
    Ganma2 = Ganma2[0:2*N:]
    Theta = np.append(Theta, Ganma2, 1)
    
Ganma = np.reshape(Ganma,(2*N,1))

#二次計画問題
Q = np.zeros((2*N, 2*N))
R = np.zeros((N, N))
for i in range(2*N):
    Q[i][i] = 1
for i in range(N):
    R[i][i] = 1
    
#Pのこと
H = np.dot(np.dot(Theta.T, Q), Theta) + R
#gbtなに？？
gbt = -2*np.dot(Theta.T, Q)

#拘束条件
umax = 1
umin = -1
F = np.zeros((2*N,N+1))
for i in range(N):
    F[2*i][i]=1
    F[2*i+1][i]=-1
    F[2*i][N]=-1*umax
    F[2*i+1][N]=umin
    
FF=np.zeros((2*N,N+1))
for i in range(2*N):       #制約の数
    for j in range(N):
        hoge = 0
        for k in range(j,N):
            hoge = hoge + F[i][k]

        FF[i][j] = hoge

    FF[i][N] = F[i][N]    #定数部
FFu = FF[:,0:N]
FFp = -1*FF[:,N]
FFp = FFp.reshape(2*N,1)

#シュミレーション開始
#x(k)
xk_plant = np.array([[1], [0]])
#u(k-1)
uk_1 = 0
#最終目標値
set_point = np.array([[0], [0]])

save_x = np.array(xk_plant)
#なんでいる？？
save_u = np.array([uk_1])

Time = 100
for i in range(Time):

    """
    目標値ベクトル  T
    """
    lamda  = 0.9
    epsilonk_1 = lamda * (set_point - xk_plant)
    T=np.zeros((2*N,1))
    for i in range(N):
        T[2*i,0]= (set_point - epsilonk_1)[0,0]
        T[2*i+1,0]= (set_point - epsilonk_1)[1,0]
        epsilonk_1 = lamda * epsilonk_1


    """
    誤差ベクトル
    """
    Epsilon = T - np.dot(Phi,xk_plant) - np.dot(Ganma,uk_1)
    """
    ソルバーにより解く
    入力決定
    """
    g=np.dot(gbt,Epsilon)
    gm = g
    Hm = H
    FFum = FFu
    FFpm = FFp - uk_1*FF[:,0].reshape(2*N,1)

    x = cvxpy.Variable(N)
    objective = cvxpy.Minimize(0.5 * cvxpy.quad_form(x, Hm) + gm.T @ x)
    #constraints = [np.dot(FFum, xk_plant) <= FFpm]
    problem = cvxpy.Problem(objective)
    sol = problem.solve()

    uk = sol+uk_1

    """
    モデルに入力を加える
    """

    xk_1_plant = np.dot(Ad,xk_plant) + np.dot(Bd,uk)

    """
    保存
    """
    save_x = np.append(save_x,xk_1_plant,1)
    save_u = np.append(save_u,np.array([uk]),0)

    """
    更新
    """
    xk_plant = xk_1_plant
    uk_1 = uk





#グラフ化
fig = plt.figure(figsize=(12, 8)) #...1

# Figure内にAxesを追加()
ax = fig.add_subplot(211) #...2
ax2 = fig.add_subplot(212) #...2
t=np.arange(Time+1)
ax.plot(t,save_x[0,:])
ax.plot(t,save_x[1,:])
ax2.plot(save_u)
plt.show()


