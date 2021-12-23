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
C=np.array([1,0])
D=np.array([0])

dt = 0.01
dsys = scipy.signal.cont2discrete((A,B,C,D),dt,method ='zoh')

Ad = dsys[0]
Bd = dsys[1]
Cd = dsys[2]

print("----dicred-----")
print(Ad)
print(Bd)
print(Cd)


umax = 1
umin = -1

"""
レギュレータ　リカッチを解く
"""
q=1
Qq = np.array([[q,0],[0 ,q]])
Rr = np.array([1])

X = scipy.linalg.solve_discrete_are(Ad,Bd,Qq,Rr)
print(X)
hoge = np.dot(np.dot(Ad.T,X),Ad) + Qq 
hoge2 = np.dot(np.dot(Ad.T,X),Bd)
hoge3 = np.linalg.inv(Rr+np.dot(np.dot(Bd.T,X),Bd))
hoge = hoge +  np.dot(np.dot ( hoge2 ,hoge3 ) ,hoge2.T) 
print(hoge)

K=-1*np.dot(hoge3,hoge2.T)
print("---K----")
print(K)

#シミュレーション開始
xk_plant = np.array([[1],[0]])  #x(k)
uk_1 = 0                        #u(k-1)
set_point =  np.array([[0],[0]]) #最終目標値

save_x = np.array(xk_plant)
save_u = np.array([uk_1])
#save_r


Time = 1000                    #Time＊dt  [s]
for i in range(Time):

    uk = np.dot(K,xk_plant)[0,0]

    if uk>umax:
        uk = umax
    elif uk<umin:
        uk = umin



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