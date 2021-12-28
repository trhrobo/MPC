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
//#include "matplotlibcpp.h"

//namespace plt = matplotlibcpp;

/*
非線形最適制御入門(例8.1)
*/

Eigen::Matrix<double, 2, 1> calcModelF(Eigen::Matrix<double, 2, 1> _X, double _u, double _t){
    Eigen::Matrix<double, 2, 1> model_F;
    double x1 = _X(0, 0);
    double x2 = _X(1, 0);
    model_F << x2,
         (1-std::pow(x1, 2)-std::pow(x2, 2))*x2-x1+_u;
    return model_F;
}

Eigen::Matrix<double, 2, 1> caclRphiRx(Eigen::Matrix<double, 2, 1> _X, double _u, Eigen::Matrix<double, 2, 1> _lamda, double _t){
    Eigen::Matrix<double, 2, 1> RphiRx;
    double _x1_f = _X(0, 0);
    double _x2_f = _X(1, 0);
    RphiRx << _x1_f,
              _x2_f;
    return RphiRx;
}

Eigen::Matrix<double, 30, 1> GMRES(double A, double b){
    //Ax = bの連立一次方程式xについてを解く
    Eigen::Matrix<double, 30, 1> dU;

}

int main(){
    /*各種定数設定*/
    //目標値に対する誤差
    constexpr double error=0.00001;
    constexpr double dt=0.001;
    //予測ステップ
    constexpr int N_step=30;
    constexpr double T_predict=10;
    constexpr double dtau=T_predict/N_step;
    constexpr double zeta=0.1;
    constexpr double h=0.001;
    //初期値設定
    double t=0;
    //各時刻における制御入力
    double u=0;
    constexpr double x1=10;
    constexpr double x2=10;
    //x_(x*)
    Eigen::Matrix<double, N_step, 2> X_=Eigen::MatrixXd::Ones(N_step, 2);
    Eigen::Matrix<double, N_step, 2> Lamda_=Eigen::MatrixXd::Ones(N_step, 2);
    //1
    //X(0)を測定する
    Eigen::Matrix<double, 2, 1> X;
    X << x1,
         x2;
    //U(0)を決定する
    Eigen::Matrix<double, N_step, 1> U=Eigen::MatrixXd::Ones(N_step, 1);
    while(1){
        //2
        //u(t)=u0(t)をシステムへの制御入力とする
        u=U(0, 0);
        //3
        //x_を求める
        //x0*(t)=x(t)
        X_(0, 0)=X(0, 0);
        X_(0, 1)=X(1, 0);
        Eigen::Matrix<double, 2, 1> prev_temp_X_;
        prev_temp_X_ << X_(0, 0),
                        X_(0, 1);
        for(int i=1; i <= (N_step-1); i++){
            Eigen::Matrix<double, 2, 1> temp_X_;
            temp_X_=prev_temp_X_+calcModelF(prev_temp_X_, U(i, 0), t+i*dtau)*dtau;
            X_(i, 0)=temp_X_(0, 0);
            X_(i, 1)=temp_X_(1, 0);
            prev_temp_X_=temp_X_;
        }
        //4
        //lamda_を求める
        //lamdaN_を決定する
        Eigen::Matrix<double, 2, 1> LamdaN_=caclRphiRx(X_(N_step-1, 0), X_(N_step-1, 1));
        Lamda_(N-1, 0)=LamdaN_(0, 0);
        Lamda_(N-1, 1)=LamdaN_(1, 0);
        //lamda_を計算する
        Eigen::Matrix<double, 2, 1> prev_temp_Lamda_;
        prev_temp_Lamda_ << Lamda_(N-1, 0),
                            Lamda_(N-1, 1);
        for(int i=N-1; i > 0; --i){
            Eigen::Matrix<double, 2, 1> temp_Lamda_;
            //FIXME:X_の行列を転置行列にすればいい
            temp_Lamda_=prev_temp_Lamda_+caclRphiRx(X_[i], U[i], prev_temp_Lamda_, t+i*dtau)*dtau;
            Lamda_(i, 0)=temp_Lamda_(0, 0);
            Lamda_(i, 1)=temp_Lamda_(1, 0);
        }
        //5
        //gmres法を用いてdUを求める
        //Fを作る
        F=createF();
        dU(t)=GMRES();
        U(t+dt)=U(t)+dU(t)*dt;
        /*----------------------------------------------------------------
        ------------------------------------------------------------------*/
        //Ax=bを解く
        constexpr double m=30;
        Eigen::Matrix<double, 30, 1> gmres_X;
        std::array<double, 30> gmres_V;
        Eigen::Matrix<double, 30, 1> gmres_G;
        Eigen::Matrix<double, 30, 1> gmres_X0=Eigen::MatrixXd::Zero(30, 1);
        //初期残差を測定
        Eigen::Matrix<double, 30, 1> gmres_R0=Eigen::MatrixXd::Zero(30, 1);
        //gmres_r0を正規化
        //1
        //FIXME:行列に割り算はない
        gmres_X[0]=gmres_X0/gmres_X0.norm();
        //2
        //初期残差の測定
        gmres_R0=b-A*gmres_X0;
        //FIXME:行列に割り算はない
        gmres_V[1]=gmres_R0/gmres_R0.norm();
        //3
        gmres_G[0]=gmres_R0.norm();
        //4
        for(int j=1; j<=m; ++j){
            //5
            //6
            //7
            //8
            //9
            //10
            for(int i=0; i<=(j-1); ++i){
                //11
                //12
                //13
                //14
                //15
            }
            //16
            //17
            //18
            //19
            //20
            //21
        }
        //22
        //23
        //24
        //25
        //26
        //27
        //28
        /*----------------------------------------------------------------
        ------------------------------------------------------------------*/
        //6
        //x(t)=x(t+dt)でxの更新
        //FIXME:calcdXの引数に直接temp_x1, temp_x2を入れたら綺麗になる
        double temp_x1=X(0, 0);
        double temp_x2=X(0, 1);
        X+=calcModelF(temp_x1, temp_x2, u)*dt;

        //7
        t=t+dt
        if((goal-val) < error){
            break;
        }
    }
}