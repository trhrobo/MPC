#include<iostream>
#include<cmath>
#include<vector>
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
constexpr double error=0.001;
//予測ステップ
constexpr int N_step=10;
constexpr double alpha=0.5;
constexpr double tf=1.0;
constexpr double zeta=100.0;
constexpr double h=0.01;
//初期値設定
constexpr double x1=2;
constexpr double x2=0;
//mはリスタートパラメータ
constexpr int m=30;

//x={x1, x2}    
constexpr int x_size=2;
//各時刻における制御入力
//・入力が一つ
//・ダミー入力が一つ
//・等式制約条件が一つ
//u={u, v, rho}
constexpr int u_size=3;
//f={rHru(x_size), c}
constexpr int f_size=3;

class NMPC{
    public:
        NMPC(Eigen::Matrix<double, x_size, 1> _X){
            X=_X;
            for(int i=0; i<N_step; ++i){
                U(u_size*i, 0)=0;
                U((u_size*i)+1, 0)=0.49;
                U((u_size*i)+2, 0)=0.011;
            }
        }
        void updateState(Eigen::Matrix<double, u_size, 1> _u){
            //FIXME:ルンゲクッタを解く
            //状態Xを更新する
            X=;
        }
        Eigen::Matrix<double, x_size, 1> calModel(Eigen::Matrix<double, x_size, 1> _X, Eigen::Matrix<double, u_size, 1> _U, double _t){
            Eigen::Matrix<double, x_size, 1> model_F;
            double x1=_X(0, 0);
            double x2=_X(1, 0);
            //U={u, v, rho}
            double u=_U(0, 0);
            model_F << x2,
                       ((1-x1*x1-x2*x2)*x2-x1+u);
            return model_F;
        }
        Eigen::Matrix<double, 1, x_size> rphirx(Eigen::Matrix<double, x_size, 1> _X, double _t){
            Eigen::Matrix<double, 1, x_size> rphirx=_X.transpose();
            return rphirx;
        }
        Eigen::Matrix<double, 1, x_size> rHru(Eigen::Matrix<double, x_size, 1> _x_, Eigen::Matrix<double, u_size, 1> _u_, Eigen::Matrix<double, x_size, 1> _lamda_){
            double u=_u_(0, 0);
            double rho=_u_(1, 0);
            double lamda2=_lamda_(1, 0); 
            //vはダミー変数
            double v=std::sqrt(0.5*0.5-u*u);
            Eigen::Matrix<double, 1, x_size> ans;
            ans<<u+lamda2+2*rho*u, -0.01+2*rho*v;
            return ans;
        }
        Eigen::Matrix<double, 1, x_size> rHrx(Eigen::Matrix<double, x_size, 1> _x_,Eigen::Matrix<double, u_size, 1> _u, Eigen::Matrix<double, x_size, 1> _lamda_, double _t){
            double x1=_x_(0, 0);
            double x2=_x_(1, 0);
            double lamda1=_lamda_(0, 0);
            double lamda2=_lamda_(1, 0);
            Eigen::Matrix<double, 1, x_size> ans;
            ans<<x1-2*x1*x2*lamda2-lamda2, x2+lamda1+(-3*x2*x2-x1*x1+1)*lamda2;
            return ans;
        }
        Eigen::Matrix<double, x_size, 1> Constraint(double _u_, double _x_, double _lamda_, double _rho_){
            //制約なし
            Eigen::Matrix<double, x_size, 1> ans=Eigen::MatrixXd::Zero(x_size, 1);
            return ans;
        }
        Eigen::Matrix<double, u_size, 1> CGMRES(double _time){
            dt=tf*(1-std::exp(-alpha*_time))/N_step;
            //gmres法を用いてdUを求める
            Eigen::Matrix<double, u_size*N_step, 1> dU=Eigen::MatrixXd::Zero(u_size*N_step, 1);
            Eigen::Matrix<double, u_size*N_step, 1> gmres_Xm=Eigen::MatrixXd::Zero(u_size*N_step, 1);
            Eigen::Matrix<double, u_size*N_step, 1> gmres_X0=Eigen::MatrixXd::Zero(u_size*N_step, 1);
            Eigen::Matrix<double, m, 1> gmres_Ym=Eigen::MatrixXd::Zero(m, 1);
            Eigen::Matrix<double, u_size*N_step, 1> gmres_V[m];
            Eigen::Matrix<double, u_size*N_step, m> gmres_Vm=Eigen::MatrixXd::Zero(u_size*N_step, m);
            Eigen::Matrix<double, u_size*N_step, 1> gmres_R0=Eigen::MatrixXd::Zero(u_size*N_step, 1);
            //gmres_X0を設定する
            gmres_X0=U;
            //初期残差gmres_R0を求める
            gmres_R0=calR0(U, X, t);
            gmres_V[0]=gmres_R0.normalized();
            double g[m+1]{};
            g[0]=gmres_R0.norm();
            double h[m+1][m]{};
            double r[m][m]{};
            double c[m]{};
            double s[m]{}; 
            for(int i=0; i<m; ++i){
                for(int k=0; k<(i+1); ++k){
                    Eigen::Matrix<double, u_size*N_step, 1> tempAv=calAv(U, X, gmres_V[i], t);
                    h[k][i]=tempAv.dot(gmres_V[k]);
                }
                Eigen::Matrix<double, u_size*N_step, 1> temp_sigma=Eigen::MatrixXd::Zero(u_size*N_step, 1);
                for(int k=0; k<(i+1); k++){
                    temp_sigma+=h[k][i]*gmres_V[k];
                }
                Eigen::Matrix<double, u_size*N_step, 1>temp_V=calAv(U, X, gmres_V[i], t)-temp_sigma;
                h[i+1][i]=temp_V.norm();
                gmres_V[i+1]=temp_V/h[i+1][i];
                r[0][i]=h[0][i];
                for(int k=0; k<((i+1)-1); ++k){
                    double temp1=c[k]*r[k][i]+s[k]*h[k+1][i];
                    double temp2=-1*s[k]*r[k][i]+c[k]*h[k+1][i];
                    r[k][i]=temp1;
                    r[k+1][i]=temp2;
                }
                c[i]=r[i][i]/std::sqrt(r[i][i]*r[i][i]+h[i+1][i]*h[i+1][i]);
                s[i]=h[i+1][i]/std::sqrt(r[i][i]*r[i][i]+h[i+1][i]*h[i+1][i]);
                g[i+1]=-1*s[i]*g[i];
                g[i]=c[i]*g[i];
                r[i][i]=c[i]*r[i][i]+s[i]*h[i+1][i];
            }
            Eigen::Matrix<double, m, 1> Gm;
            for(int i=0; i<m; ++i){
                Gm(i, 0)=g[i];
            }
            Eigen::Matrix<double, m, m> Rm;
            for(int i=0; i<m; ++i){
                for(int k=0; k<m; ++k){
                    Rm(i, k)=r[i][k];
                }
            }
            //上で作ったgmres_V[i]をgmres_Vmに代入する
            for(int i=0; i<m; ++i){
                gmres_Vm.col(i)=gmres_V[i];
            }
            //後退代入によってRm*Ym=Gmを解く
            for(int i=(m-1); i>=0; --i){
                double temp_sigma_back=0;
                for(int k=(i+1); k<m; ++k){
                    temp_sigma_back+=Rm(i, k)*gmres_Ym[k];
                }
                gmres_Ym[i]=(Gm[i]-temp_sigma_back)/Rm(i, i);
            }
            gmres_Xm=gmres_X0+gmres_Vm*gmres_Ym;
            dU=gmres_Xm;
            U+=dU*dt;
            return U.block(0, 0, u_size, 1);
        }
        Eigen::Matrix<double, x_size*N_step, 1> calF(Eigen::Matrix<double, x_size, 1> _x, double _t){
            Eigen::Matrix<double, x_size*N_step, 1> F;
            //制約なし
            //0~Nまでx1, x2, lamda1, lamda2
            Eigen::Matrix<double, x_size*N_step, 1> X_;
            Eigen::Matrix<double, x_size*N_step, 1> Lamda_;
            //x_(x*)を求める
            //x0*(t)=x(t)を計算する
            X_.block(0, 0, x_size, 1)=_x;
            //xi+1*=xi*+f(xi*,ui*,t+i*dtau)*dtauを計算する
            Eigen::Matrix<double, x_size, 1> prev_X_=_x;
            for(int i=1; i < N_step; i++){
                X_.block(x_size*i, 0, x_size, 1)=prev_X_+calModel(prev_X_, U.block(u_size*i, 0, u_size, 1), _t+i*dt)*dt;
                prev_X_=X_.block(x_size*i, 0, x_size, 1);
            }
            //4
            //lamda_(lamda*)を求める
            //lamdaN*=(rphi/rx)^T(xN*,t+T)を計算する
            Eigen::Matrix<double, 1, x_size> temp_rphirx=rphirx(X_.block(x_size*(N_step-1), 0, x_size, 1), _t+T_predict);
            Lamda_.block(x_size*(N_step-1), 0, x_size, 1)=temp_rphirx.transpose();
            //lamdai*=lamdai+1*+(rH/ru)^T*dtau
            Eigen::Matrix<double, x_size, 1> prev_Lamda_=Lamda_.block(x_size*(N_step-1), 0, x_size, 1);
            //逆順で解く
            //N_step-2の理由(N_step-1で最後のLamdaのグループなので(上でそこは計算してる),それの前だからN-2)
            for(int i=(N_step-2); i >= 0; --i){
                Eigen::Matrix<double, 1, x_size> temp_rHrx=rHrx(X_.block(x_size*i, 0, x_size, 1), U.block(u_size*i, 0, u_size, 1), prev_Lamda_, _t+i*dt);
                Lamda_.block(x_size*i, 0, x_size, 1)=prev_Lamda_+temp_rHrx.transpose()*dtau;
                prev_Lamda_=Lamda_.block(x_size*i, 0, x_size, 1);
            }
            for(int i=0; i<N_step; i++){
                double lam_2=Lamda_((i*x_size)+1, 1);
                double u=U(i*u_size, 0);
                double v=U((i*u_size+1), 0);
                double rho=U((i*u_size+2), 0);
                F((i*f_size),   0)=u+lam_2+2.0*rho*u;
                F((i*f_size)+1, 0)=-0.01+2.0*rho*v;
                F((i*f_size)+2, 0)=u*u+v*v-0.5*0.5;
            }
            return F;
        }   
        Eigen::Matrix<double, x_size*N_step, 1> calAv(Eigen::Matrix<double, x_size, 1> _X, Eigen::Matrix<double, x_size*N_step, 1> _V, double _t){
            Eigen::Matrix<double, u_size, 1> tempU= _U.block(0, 0, u_size,1);
            Eigen::Matrix<double, x_size*N_step, 1> temp1=calF(U+h*_V, _X+h*calModel(_X, tempU, 0), t+h);
            Eigen::Matrix<double, x_size*N_step, 1> temp2=calF(U,      _X+h*calModel(_X, tempU, 0), t+h);
            Eigen::Matrix<double, x_size*N_step, 1> Av=(temp1-temp2)/h;
            return Av;
        }
        Eigen::Matrix<double, x_size*N_step, 1> calR0(Eigen::Matrix<double, x_size, 1> _X, double _t){
            //U'(0)=U0を使用する
            Eigen::Matrix<double, u_size*N_step, 1> dU=U;
            Eigen::Matrix<double, u_size, 1> tempU= U.block(0, 0, u_size,1);
            Eigen::Matrix<double, x_size*N_step, 1> temp11=calF(U, _X+h*calModel(_X, tempU, 0), t+h);
            Eigen::Matrix<double, x_size*N_step, 1> temp12=calF(U, _X, t);
            Eigen::Matrix<double, x_size*N_step, 1> temp1=(temp11-temp12)/h;
            Eigen::Matrix<double, x_size*N_step, 1> temp21=calF(U+h*dU, _X+h*calModel(_X, tempU, 0), t+h);
            Eigen::Matrix<double, x_size*N_step, 1> temp22=calF(U, _X+h*calModel(_X, tempU, 0), t+h);
            Eigen::Matrix<double, x_size*N_step, 1> temp2=(temp21-temp22)/h;
            Eigen::Matrix<double, x_size*N_step, 1> R0=-1*zeta*calF(U, _X, 0)-temp1-temp2;
            return R0;
        }
    private:
        Eigen::Matrix<double, x_size, 1> X=Eigen::MatrixXd::Zero(x_size, 1);
        Eigen::Matrix<double, u_size*N_step, 1> U=Eigen::MatrixXd::Zero(u_size*N_step, 1);
        double dt;
};
int main(){
    constexpr double dt=0.01;
    constexpr int iteration_time=20;
    constexpr int iteration_num=iteration_time/dt;
    //現在の状態
    //x={x1, x2}
    Eigen::Matrix<double, 2, 1> initX;
    //X(0)を測定する(初期値を代入する)
    initX << x1,
             x2;
    NMPC nmpc(initX);
    for(int i=0; i<iteration_num; ++i){
        double time=i*dt;
        Eigen::Matrix<double, u_size, 1> u=nmpc.CGMRES(time);
        nmpc.updateState(u);
    }
}