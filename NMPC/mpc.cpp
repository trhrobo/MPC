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

/*各種定数設定*/
//目標値に対する誤差
constexpr double error=0.00001;
constexpr double dt=0.001;
//予測ステップ
constexpr int N_step=10;
constexpr double T_predict=10;
constexpr double dtau=T_predict/N_step;
constexpr double zeta=0.1;
constexpr double h=0.001;
//初期値設定
double t=0;
constexpr double x1=10;
constexpr double x2=10;
//mはリスタートパラメータ
constexpr int m=30;
Eigen::Matrix<double, 2, 1> calModel(Eigen::Matrix<double, 2, 1> _X, Eigen::Matrix<double, 2, 1> _U, double _t){
    Eigen::Matrix<double, 2, 1> model_F;
    double x1=_X(0, 0);
    double x2=_X(1, 0);
    //U={u, rho}
    double u=_U(0, 0);
    model_F << x2,
               ((1-x1*x1-x2*x2)*x2-x1+u);
    return model_F;
}
Eigen::Matrix<double, 2, 1> rphirx(Eigen::Matrix<double, 2, 1> _X, double _t){
    Eigen::Matrix<double, 2, 1> rphirx=_X;
    return rphirx;
}
Eigen::Matrix<double, 2, 1> rHru(Eigen::Matrix<double, 2, 1> _x_, Eigen::Matrix<double, 2, 1> _u_, Eigen::Matrix<double, 2, 1> _lamda_){
    double u=_u_(0, 0);
    double rho=_u_(1, 0);
    double lamda2=_lamda_(1, 0);
    //vはダミー変数
    double v=std::sqrt(0.5*0.5-u*u);
    Eigen::Matrix<double, 2, 1> ans;
    ans<<u+lamda2+2*rho*u
         -0.01+2*rho*v;
    return ans;
}
Eigen::Matrix<double, 2, 1> rHrx(Eigen::Matrix<double, 2, 1> _x_,Eigen::Matrix<double, 3, 1> _u, Eigen::Matrix<double, 2, 1> _lamda_, double _t){
    double x1=_x_(0, 0);
    double x2=_x_(1, 0);
    double lamda1=_lamda_(0, 0);
    double lamda2=_lamda_(1, 0);
    Eigen::Matrix<double, 2, 1> ans;
    ans<<x1-2*x1*x2*lamda2-lamda2,
         x2+lamda1+(-3*x2*x2-x1*x1+1)*lamda2;
    return ans;
}
Eigen::Matrix<double, 2, 1> Constraint(double _u_, double _x_, double _lamda_, double _rho_){
    //制約なし
    Eigen::Matrix<double, 2, 1> ans=Eigen::MatrixXd::Zero(2, 1);
    return ans;
}
Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> calF(Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _U, Eigen::Matrix<double, 2, 1> _x, double _t){
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> F;
    //制約なし
    //0~Nまでx1, x2, lamda1, lamda2
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> X_;
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> Lamda_;
    //x_(x*)を求める
    //x0*(t)=x(t)を計算する
    X_(0, 0)=_x;
    //xi+1*=xi*+f(xi*,ui*,t+i*dtau)*dtauを計算する
    Eigen::Matrix<double, 2, 1> prev_temp_X_=X_(0, 0);
    for(int i=1; i < N_step; i++){
        X_(i, 0)=prev_temp_X_+calModel(prev_temp_X_, _U(i, 0), t+i*dtau)*dtau;
        prev_temp_X_=X_(i, 0);
    }
    //4
    //lamda_(lamda*)を求める
    //lamdaN*=(rphi/rx)^T(xN*,t+T)を計算する
    Lamda_(N_step-1, 0)=rphirx(X_(N_step-1, 0), t+T_predict);
    //lamdai*=lamdai+1*+(rH/ru)^T*dtau
    Eigen::Matrix<double, 2, 1> prev_temp_Lamda_=Lamda_(N_step-1, 0);
    //逆順で解く
    for(int i=N_step-1; i > 0; --i){
        Eigen::Matrix<double, 2, 1> temp_Lamda_;
        //FIXME:X_の行列を転置行列にすればいい
        Lamda_(i, 0)=prev_temp_Lamda_+rHrx(X_(i, 0), _U(i, 0), prev_temp_Lamda_, t+i*dtau)*dtau;
        prev_temp_Lamda_=Lamda_(i, 0);
    }
    for(int i=0; i<N_step; i++){
        F(i, 0)=rHru(X_(i, 0), _U(i, 0), Lamda_(i, 0));
    }
    return F;
}
/*
Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> calAv(Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _U_, Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _X_, Eigen::Matrix<double, n, 1> _V){
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> Av=(calF(_U_+h*_V_, x+h*dx, 0)-calF(_U_, x+h*dx, 0))/h;
    return Av;
}
*/
Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> calAv(Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _U, Eigen::Matrix<double, 2, 1> _X, Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _V){
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> Av=(calF(_U+h*_V, _X+h*calModel(_X, _U(0, 0), 0))-calF(_U, _X+h*calModel(_X, _U(0, 0), 0))/h;
    return Av;
}
Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> calR0(Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _U_, Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> _X_){
    //U'(0)=U0を使用する
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> R0=-1*zeta*calF(_U_, _X_, 0) -(calF(_U_, _X_+hdx, h)-F(_U_, _X_, 0))/h-(F(_U_+h*dU, X+h*dX, h)-calF(_U_, _X_+h*dx, h))/h;
    return R0;
}
int main(){
    //現在の状態
    //x={x1, x2}
    Eigen::Matrix<double, 2, 1> X;
    //X(0)を測定する(初期値を代入する)
    X << x1,
         x2;
    //U(0)を決定する
    //各時刻における制御入力
    //・入力が一つ
    //・ダミー入力が一つ
    //・等式制約条件が一つ
    constexpr int u_size=3;
    //u={u, dummy, rho}
    //dummyは入れなくていいのでは?
    Eigen::Matrix<double, 3, 1> u=Eigen::MatrixXd::Zero(3, 1);
    //FIXME:Uの初期化
    Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, 1> U;
    while(1){
        //u(t)=u0(t)をシステムへの制御入力とする
        //u={u, dummy, rho}
        u=U(0, 0);
        //gmres法を用いてdUを求める
        Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, 1> dU;
        /*----------------------------------------------------------------
        ------------------------------------------------------------------*/
        Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, 1> gmres_Xm;
        Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, 1> gmres_X0;
        Eigen::Matrix<double, m, 1> gmres_Ym;
        Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, 1> gmres_V[m];
        Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, m> gmres_Vm;
        Eigen::Matrix<Eigen::Matrix<double, 3, 1>, N_step, 1> gmres_R0;
        //初期残差gmres_R0を求める
        gmres_R0=calR0();
        gmres_V[0]=gmres_R0.normalized();
        //Vmを作る
        double h[m][m]{};
        for(int i=0; i<m; ++i){
            for(int k=0; k<m; ++k){
                h[k][i]=calAv(gmres_V[i]).dot(gmres_V[k]);
            }
            Eigen::Matrix<double, N_step, 1> temp_sigma=Eigen::MatrixXd::Zero(N_step, 1);
            for(int k=0; k<i; k++){
                temp_sigma=h[i][k]*gmres_V[i];
            }
            Eigen::Matrix<double, n, 1>temp_V=calAv(gmres_V[i])-temp_sigma;
            double temp_size_V=temp_V.norm();
            gmres_V[i]=(1.0/temp_size_V)*temp_V;
        }
        for(int i=0; i<n; ++i){
            for(int k=0; k<m; ++k){
                //FIXME:代入を間違えているかも
                gmres_Vm(i, k)=gmres_V[i][k];
            }
        }
        //最小化問題を解く
        //FIXME:c, s, rを求める必要がある
        double c[m]{};
        double s[m]{}; 
        double r[m]{};
        //Rmを作る
        Eigen::Matrix<double, m, m> Hm;
        Eigen::Matrix<double, m+1, m> _Hm;
        for(int i=0; i<m; ++i){
            for(int k=0; k<m; ++k){
                //FIXME:代入を間違えているかも
                Hm(i, k)=h[i][k];
            }
        }
        Eigen::Matrix<double, 1, m> temp_Hm=Eigen::MatrixXd::Zero(1, m);
        //FIXME:代入法を間違えているかも
        temp_Hm(0, m-1)=h[m][m-1];
        _Hm = (Hm,
               temp_Hm);
        //Givens回転を用いて_Hmを上三角行列に変換する
        Eigen::Matrix<double, m+1, m+1> Qm=Eigen::MatrixXd::Identity(m+1, m+1);
        for(int i=0; i<m; i++){
            Eigen::Matrix<double, m+1, m+1> Omega=Eigen::MatrixXd::Identity(m+1, m+1);
            //TODO:回転行列をその部分に入れるようにした方が綺麗
            Omega(i, i)=c[i];
            Omega(i+1, i)=-1*s[i];
            Omega(i, i+1)=s[i];
            Omega(i+1, i+1)=c[i];
            Qm*=Omega;
        }
        Eigen::Matrix<double, m+1, m> _Rm;
        _Rm=Qm*_Hm;
        //FIXME:_Rmの最下段を取り除いたものをRmとする必要がある
        Eigen::Matrix<double, m, m> Rm=_Rm;
        //gmを作る
        double g[m]{};
        Eigen::Matrix<double, m, 1> Gm;
        for(int i=0; i<m; i++){
            double temp_prod_s=1;
            for(int k=0; k<(i-1); k++){
                temp_prod_s*=s[k];
            }
            //FIXME:std::pow(-1, i-1)を偶数,奇数で判別した方が速くなる
            g[i]=std::pow(-1, i-1)*gmres_R0.norm()*c[i]*temp_prod_s;
            Gm(0, i)=g[i];
        }
        //後退代入によってRm*Ym=Gmを解く
        double temp_sigma_back;
        for(int i=m; i>0; --i){
            double temp_sigma_back=0;
            for(int k=0; k<(m-k); ++k){
                temp_sigma_back+=Gm[k]*gmres_Ym[k];
            }
            gmres_Ym[i]=(Gm[i]-temp_sigma_back)/Rm(i, i);
        }
        gmres_Xm=gmres_X0+gmres_Vm*gmres_Ym;
        dU=gmres_Xm;
        /*----------------------------------------------------------------
        ------------------------------------------------------------------*/
        U=U+dU*dt;
        //x(t)=x(t+dt)でxの更新
        //FIXME:calcdXの引数に直接temp_x1, temp_x2を入れたら綺麗になる
        double temp_x1=X(0, 0);
        double temp_x2=X(0, 1);
        X+=calcModel(temp_x1, temp_x2, u)*dt;
        t=t+dt
        if((goal-val) < error){
            break;
        }
    }
}