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
constexpr double zeta = 0.1;
constexpr double h = 0.001;

Eigen::Matrix caldX(double _x1, double _x2, double _u){
    Eigen::Matrix<double, 2, 1> _dX;
    _dX << _x2,
            (1-std::pow(_x1, 2)-std::pow(_x2, 2))*_x2-_x1+_u;
    return _dX;
}
void Hamilton()
void GMRES(double A, double b){
    //Ax = bの連立一次方程式xについてを解く
}

int main(){
    //モデル
    Matrix<double, 2, 1> dX = caldX();

    //初期値設定
    double t = 0;
    constexpr double x1 = 10;
    constexpr double x2 = 10;
    //1
    //X(0)を測定する
    Eigen::Matrix<double, 2, 1> X;
    X << x1,
         x2;
    //U(0)を決定する
    Eigen::Matrix<double, N_step, 1>;
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
        //Fを作る
        F = createF();
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