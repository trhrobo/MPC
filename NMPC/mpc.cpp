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
        //mはリスタートパラメータ
        constexpr int m=30;
        //FIXME:nを求める必要がある
        constexpr int n=40;
        //num_solutionは解の個数
        constexpr double num_solution=60;
        Eigen::Matrix<double, n, m> A;
        Eigen::Matrix<double, n, 1> b;
        Eigen::Matrix<double, num_solution, 1> gmres_Xm;
        Eigen::Matrix<double, num_solution, 1> gmres_X0=Eigen::MatrixXd::Zero(num_solution, 1);
        Eigen::Matrix<double, n, 1> gmres_V[m];
        Eigen::Matrix<double, n, m> gmres_Vm;
        Eigen::Matrix<double, n, 1> gmres_R0;
        gmres_R0=b-A*gmres_X0;
        gmres_V[0]=gmres_R0.normalized();
        //Vmを作る
        double h[m][m]{};
        for(int i=0; i<m; ++i){
            for(int k=0; k<m; ++k){
                //FIXME:サイズがおかしい
                //FIXME:hの求め方を間違えてるかも
                h[k][i]=A*gmres_V[i].dot(gmres_V[k]);
            }
            Eigen::Matrix<double, n, 1> temp_sigma=Eigen::MatrixXd::Zero(n, 1);
            for(int k=0; k<i; k++){
                temp_sigma=h[i][k]*gmres_V[i];
            }
            Eigen::Matrix<double, n, 1>temp_V=A*gmres_V[i]-temp_sigma;
            double temp_size_V=temp_V.norm();
            gmres_V[i]=(1.0/temp_size_V)*temp_V;
        }
        for(int i=0; i<n; ++i){
            for(int k=0; k<m; ++k){
                //FIXME:代入を間違えている
                gmres_Vm[i][k]=gmres_V[i][k];
            }
        }
        //最小化問題を解く
        double c[m]{};
        double s[m]{}; 
        for(int j=1; j<=m; ++j){
            //Arnoldi法を行う
            for(int i=1; i<=j; ++i){
                h[i][j]=A*v[j].dot(v[i]);
            }
            double temp_sigma;
            for(int i=1; i<=j; ++i){
                temp_sigma=h[i][j]*v[i];
            }
            _v[j+1]=A*v[j] - temp_signa;
            h[j+1][j]=_v[j+1].norm();
            v[j+1]=_v[j+1]/h[j+1][j];
            //FIXME:(0)はどうしたらいい?
            gmres_R[1][j]=h[1][j];
            for(int i=1; i<=(j-1); ++i){
                //FIXME:(i-1)はどうしたらいい?
                double temp1=c[i]*gmres_R[i][j]+s[i]*h[i+1][j];
                double temp2=-1*s[i]*gmres_R[i][j]+c[i]*h[i+1][j];
                gmres_R[i][j]=temp1;
                gmres_R[i+1][j]=temp2;
            }
            c[j]=gmres_R[j][j]/std::sqrt(std::pow(gmres_R[j][j], 2)+std::pow(h[j+1][j], 2));
            s[j]=h[j+1][j]/std::sqrt(std::pow(gmres_R[j][j], 2)+std::pow(h[j+1][j], 2));
            gmres_R[j][j]=c[j]*gmres_R[j][j]+s[j]*h[j+1][j];
        }
        Eigen::Matrix<double, 30> Ym=Rm.colPivHouseholderQr().solve(Gm);
        Eigen::Matrix<double, 30> Xm=X0+Vm*Ym;
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