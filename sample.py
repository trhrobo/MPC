"""
LMPC シミュレーション
"""
import numpy as np
import scipy
from scipy import signal 

import numpy
import cvxopt
from cvxopt import matrix

import matplotlib.pyplot as plt

c=1
k=1
A=np.array([[-c, k],[1,0]]) #不安定
B=np.array([[1],[0]])
C=np.array([0,1])
D=np.array([1])

dt = 0.01
dsys = scipy.signal.cont2discrete((A,B,C,D),dt,method ='zoh')

Ad = dsys[0]
Bd = dsys[1]
Cd = dsys[2]

print("----dicred-----")
print(Ad)
print(Bd)
print(Cd)

N=10

"""
各種定数

x(k+1)     
x(k+2)  =  Phi * x(k)  + Ganma * u(k-1) + Theta Δu(k ~ k+N-1)
:
x(k+N)

"""
Adi = Ad
Phi = Ad
for i in range(N-1):
    Adi = np.dot(Adi,Ad)
    Phi = np.append(Phi,Adi)


Phi = np.reshape(Phi,(2*N,2))
print("----Phi-----")
print(Phi)

Ganma = Bd
Adi = Ad
for i in range(N-1):
    Ganma = np.append(Ganma,np.dot(Adi,Bd))
    Adi = np.dot(Adi,Ad)

Ganma = np.reshape(Ganma,(2*N,1))
print("----Ganma-----")
print(Ganma)

zero = np.zeros((2,1))
Theta = Ganma
Ganma2 = Ganma                      #Theta作る用
for i in range(N-1):
    Ganma2 = np.append(zero,Ganma2,0)   #上にゼロをつける
    Ganma2 = Ganma2[0:2*N:]             #下を削除
    Theta = np.append(Theta,Ganma2,1)

Ganma = np.reshape(Ganma,(2*N,1))       #不必要かも
print("----Ganma-----")
print(Ganma)

"""
二次計画問題係数
"""
Q = np.zeros((2*N,2*N))
R = np.zeros((N,N))
for i in range(2*N):
    Q[i][i] = 1
for i in range(N):
    R[i][i] = 1

H=np.dot(np.dot(Theta.T,Q),Theta) + R
gbt = -2*np.dot(Theta.T,Q)

"""
拘束条件
umin < u < umax とする

F [u   < 0
   1]

FFu*u < FFp
"""
umax = 1
umin = -1
F=np.zeros((2*N,N+1))
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

#シミュレーション開始
xk_plant = np.array([[1],[0]])  #x(k)
uk_1 = 0                        #u(k-1)
set_point =  np.array([[0],[0]]) #最終目標値

save_x = np.array(xk_plant)
save_u = np.array([uk_1])
#save_r

cvxopt.solvers.options['show_progress'] = False  #計算状況　非表示
Time = 1000                    #Time＊dt  [s]
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
    gm = matrix(g)
    Hm = matrix(H)
    FFum = matrix(FFu)
    #FFpm = matrix(FFp)
    FFpm = matrix(FFp - uk_1*FF[:,0].reshape(2*N,1) )  #こっちが正しい

    sol=cvxopt.solvers.qp(Hm,gm,FFum,FFpm)

    uk = sol["x"][0]+uk_1

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