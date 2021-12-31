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
constexpr int N_step=30;
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
//FIXME:nを求める必要がある
constexpr int n=40;
Eigen::Matrix<double, 2, 1> calcModelF(Eigen::Matrix<double, 2, 1> _X, Eigen::Matrix<double, 2, 1>_u, double _t){
    Eigen::Matrix<double, 2, 1> model_F;
    double x1=_X(0, 0);
    double x2=_X(1, 0);
    double u=_X(0, 0);
    model_F << x2,
               ((1-std::pow(x1, 2)-std::pow(x2, 2))*x2-x1 _u);
    return model_F;
}

Eigen::Matrix<double, 2, 1> rphirx(Eigen::Matrix<double, 2, 1> _X, double _u, Eigen::Matrix<double, 2, 1> _lamda, double _t){
    Eigen::Matrix<double, 2, 1> rphirx;
    double _x1_f = _X(0, 0);
    double _x2_f = _X(1, 0);
    rphirx << _x1_f,
              _x2_f;
    return rphirx;
}
Eigen::Matrix<double, 2, 1> dHdu(Eigen::Matrix<double, 2, 1> _x_, Eigen::Matrix<double, 2, 1> _u_, Eigen::Matrix<double, 2, 1> _lamda_, Eigen::Matrix<double, 2, 1> _rho_){
    Eigen::Matrix<double, 2, 1> temp_lamda_=_u_;
    temp_lamda_(0, 0)=0;
    return _u_+temp_lamda_+2*_rho_*_u_;
}
Eigen::Matrix<double, 2, 1> Constraint(double _x_, double _u_, double _lamda_, double _rho_){
    //制約なし
    Eigen::Matrix<double, 2, 1> ans=Eigen::MatrixXd::Zero(2, 1);
    return ans;
}
Eigen::Matrix<double, N_step, 1> calF(Eigen::Matrix<double, 2, 1> _x_, Eigen::Matrix<double, 2, 1> _u_, Eigen::Matrix<double, 2, 1> _lamda_, Eigen::Matrix<double, 2, 1> _rho_){
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> F;
    for(int i=0; i<N_step; i++){
        F(i, 0)=dHdu(_x_, _u_, _lamda_, _rho_);
    }
    return F;
}
Eigen::Matrix<double, n, 1> calAv(Eigen::Matrix<double, n, 1> _V){
    Eigen::Matrix<double, n, 1> ans=(calF(U+h*_V, x+hx')-calF(U, x+hx'))/h;
    return ans;
}
Eigen::Matrix<double, n, 1> calA(){
    Eigen::Matrix<double, n, 1> ans;
    return ans;
}
Eigen::Matrix<double, n, m> calb(Eigen::Matrix<double, N_step, 1> _U, double _x. double _t){
    Eigen::Matrix<double, n, m> ans=-1*zeta*calF()-(calF()-calF())/h;
    return ans;
}
Eigen::Matrix<double, n, 1> calR0(){
    Eigen::Matrix<double, n, 1> ans=calb()-calA();
    return ans;
}
int main(){
    //x_(x*)
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 2> X_=Eigen::MatrixXd::Ones(N_step, 2);
    Eigen::Matrix<double, N_step, 2> Lamda_=Eigen::MatrixXd::Ones(N_step, 2);
    //1
    //X(0)を測定する
    Eigen::Matrix<double, 2, 1> X;
    X << x1,
         x2;
    //U(0)を決定する
    //各時刻における制御入力
    Eigen::Matrix<double, 2, 1> u=Eigen::MatrixXd::Zero(2, 1);
    Eigen::Matrix<Eigen::Matrix<double, 2, 1>, N_step, 1> U=Eigen::MatrixXd::Ones(N_step, 1);
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
        Eigen::Matrix<double, 2, 1> LamdaN_=rphirx(X_(N_step-1, 0), X_(N_step-1, 1));
        Lamda_(N-1, 0)=LamdaN_(0, 0);
        Lamda_(N-1, 1)=LamdaN_(1, 0);
        //lamda_を計算する
        Eigen::Matrix<double, 2, 1> prev_temp_Lamda_;
        prev_temp_Lamda_ << Lamda_(N-1, 0),
                            Lamda_(N-1, 1);
        for(int i=N-1; i > 0; --i){
            Eigen::Matrix<double, 2, 1> temp_Lamda_;
            //FIXME:X_の行列を転置行列にすればいい
            temp_Lamda_=prev_temp_Lamda_+rphirx(X_[i], U[i], prev_temp_Lamda_, t+i*dtau)*dtau;
            Lamda_(i, 0)=temp_Lamda_(0, 0);
            Lamda_(i, 1)=temp_Lamda_(1, 0);
        }
        //5
        //gmres法を用いてdUを求める
        Eigen::Matrix<double, N_step, 1> dU=Eigen::MatrixXd::Ones(N_step, 1);
        //Fを作る
        Eigen::Matrix<double, N_step, 1> F;
        for(int i=0; i<N_step; i++){
            //FIXME:サイズが違う
            F(i, 0)=dHdu();
            //制約条件無し
            //F(i+1, 0)=Constraint();
        }
        /*----------------------------------------------------------------
        ------------------------------------------------------------------*/
        //Ax=bを解く
        //num_solutionは解の個数
        constexpr double num_solution=60;
        //FIXME:サイズを間違えている可能性がある
        Eigen::Matrix<double, num_solution, 1> gmres_Xm;
        Eigen::Matrix<double, m, 1> gmres_Ym;
        Eigen::Matrix<double, num_solution, 1> gmres_X0=Eigen::MatrixXd::Zero(num_solution, 1);
        Eigen::Matrix<double, n, 1> gmres_V[m];
        Eigen::Matrix<double, n, m> gmres_Vm;
        Eigen::Matrix<double, n, 1> gmres_R0;
        //初期残差gmres_R0を求める
        gmres_R0=calR0();
        gmres_V[0]=gmres_R0.normalized();
        //Vmを作る
        double h[m][m]{};
        for(int i=0; i<m; ++i){
            for(int k=0; k<m; ++k){
                h[k][i]=calAv(gmres_V[i]).dot(gmres_V[k]);
            }
            Eigen::Matrix<double, n, 1> temp_sigma=Eigen::MatrixXd::Zero(n, 1);
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
        /*----------------------------------------------------------------
        ------------------------------------------------------------------*/
        U=dU*dt;
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