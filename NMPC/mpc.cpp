#include<iostream>
#include<cmath>
#include<limits>
#include<iomanip>
#include<vector>
#include<chrono>
#include<fstream>
#include<iomanip>
#include "/usr/include/eigen3/Eigen/Dense"
#include "/usr/include/eigen3/Eigen/Sparse"
#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/LU"
#include "matplotlibcpp.h"

/*各種定数設定*/
//目標値に対する誤差
constexpr double error = 0.00001;
constexpr double dt = 0.001;
//予測ステップ
constexpr int N_step = 30;

int main(){
    double t = 0;
    //モデル
    //1
    while(1){
        //2
        u(t) = u_(0, t);
        //3
        //x_を求める
        x_(0, t) = x(t)
        //4
        //lamda_を求める
        //5
        //gmres法を用いてdUを求める
        dU(t) = GMRES();
        U(t + dt) = U(t) + dU(t) * dt
        //6
        x(t + dt);
        //7
        t = t + dt
        if((goal - val) < error){
            break;
        }
    }
    Animation(x);
}