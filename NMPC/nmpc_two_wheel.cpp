#include<iostream>
#include<cmath>
#include<vector>
#include<chrono>
#include "/usr/include/eigen3/Eigen/Dense"
#include "/usr/include/eigen3/Eigen/Sparse"
#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/LU"
#include "/usr/include/eigen3/Eigen/QR"
#include "matplotlibcpp.h"
using namespace std::chrono;
namespace plt = matplotlibcpp;

/*
非線形最適制御入門(例8.1)
*/

/*各種定数設定*/
//目標値に対する誤差
constexpr double threshold=0.001;
//予測ステップ
constexpr int N_step=10;
constexpr double alpha=0.5;
constexpr double tf=1.0;
constexpr double zeta=100.0;
constexpr double ht=0.01;
//初期値設定
constexpr double x1=-4.5;
constexpr double x2=1.5;
constexpr double x3=0.25;
//x={x, y, theta}    
constexpr int x_size=3;
//各時刻における制御入力
//・入力が一つ
//・ダミー入力が一つ
//・等式制約条件が一つ
//u={u1, u2, v1, v2, rho1, rho2}
constexpr int u_size=6;
//f={rHru(x_size), c}
constexpr int f_size=6;
constexpr int max_iteration=u_size*N_step;
class NMPC{
    public:
        NMPC(Eigen::Matrix<double, x_size, 1> _X){
            X=_X;
            for(int i=0; i<N_step; ++i){
                U(u_size*i, 0)=1.0;
                U((u_size*i)+1, 0)=0.1;
                U((u_size*i)+2, 0)=0.1;
                U((u_size*i)+3, 0)=2.5;
                U((u_size*i)+4, 0)=0.8;
                U((u_size*i)+5, 0)=0.8;
            }
        }
        void figGraph(std::vector<double> _time){
            plt::clf();
            plt::plot(save_x1, save_x2);
            //plt::plot(_time, save_x2);
            //plt::named_plot("u", _time, save_u);
            //plt::legend();
            plt::pause(0.01);
        }
        void updateState(Eigen::Matrix<double, u_size, 1> _U, double _dt){
            //状態Xを更新する
            /*double k0[2]{};
            double k1[2]{};
            double k2[2]{};
            double k3[2]{};

            double x1=X(0, 0);
            double x2=X(1, 0);
            double u=_U(0, 0);
            save_x1.push_back(x1);
            save_x2.push_back(x2);
            save_u.push_back(u);
            for(int i=0; i<2; ++i){
                k0[i]=_dt*Func(x1, x2, u, i);
            }
            for(int i=0; i<2; ++i){
                k1[i]=_dt*Func(x1+k0[0]/2, x2+k0[1]/2, u, i);
            }
            for(int i=0; i<2; ++i){
                k2[i]=_dt*Func(x1+k1[0]/2, x2+k1[1]/2, u, i);
            }
            for(int i=0; i<2; ++i){
                k3[i]=_dt*Func(x1+k2[0], x2+k2[1], u, i);
            }
            x1+=(k0[0]+2*k1[0]+2*k2[0]+k3[0])/6;
            x2+=(k0[1]+2*k1[1]+2*k2[1]+k3[1])/6;
            X(0, 0)=x1;
            X(1, 0)=x2;*/
            X+=calModel(X, _U.block(0, 0, u_size, 1))*_dt;
            double x1=X(0, 0);
            double x2=X(1, 0);
            double x3=X(2, 0);
            double u1=_U(0, 0);
            double u2=_U(1, 0);
            save_x1.push_back(x1);
            save_x2.push_back(x2);
            save_x3.push_back(x3);
            save_u1.push_back(u1);
            save_u2.push_back(u2);
            std::cout<<x1<<" "<<x2<<" "<<x3<<std::endl;
        }
        double Func(double _x1, double _x2, double _u, double i){
            double ans{};
            if(i==0){
                ans=Func1(_x1, _x2, _u);
            }else{
                ans=Func2(_x1, _x2, _u);
            }
            return ans;
        }
        double Func1(double _x1, double _x2, double _u){
            double x1_dot=_x2;
            return x1_dot;
        }
        double Func2(double _x1, double _x2, double _u){
            double x2_dot=((1-_x1*_x1-_x2*_x2)*_x2-_x1+_u);
            return x2_dot;
        }
        Eigen::Matrix<double, x_size, 1> calModel(Eigen::Matrix<double, x_size, 1> _X, Eigen::Matrix<double, u_size, 1> _U){
            Eigen::Matrix<double, x_size, 1> model_F;
            double x1=_X(0, 0);
            double x2=_X(1, 0);
            double x3=_X(2, 0);
            //U={u, v, rho}
            double u1=_U(0, 0);
            double u2=_U(1, 0);
            model_F << std::cos(x3)*u1,
                       std::sin(x3)*u1,
                       u2;
            return model_F;
        }
        Eigen::Matrix<double, 1, x_size> rphirx(Eigen::Matrix<double, x_size, 1> _X){
            Eigen::Matrix<double, 1, x_size> rphirx=_X.transpose();
            return rphirx;
        }
        Eigen::Matrix<double, 1, x_size> rHrx(Eigen::Matrix<double, x_size, 1> _x_,Eigen::Matrix<double, u_size, 1> _u, Eigen::Matrix<double, x_size, 1> _lamda_){
            double x1=_x_(0, 0);
            double x2=_x_(1, 0);
            double x3=_x_(2, 0);
            double u1=_u(0, 0);
            double lamda1=_lamda_(0, 0);
            double lamda2=_lamda_(1, 0);
            Eigen::Matrix<double, 1, x_size> ans;
            ans<<0, 0,-lamda1*std::sin(x3)*u1+lamda2*std::cos(x3)*u1;
            return ans;
        }
        Eigen::Matrix<double, u_size, 1> CGMRES(double _time){
            dt=tf*(1-std::exp(-alpha*_time))/N_step;
            //gmres法を用いてdUを求める
            Eigen::Matrix<double, max_iteration, 1> dU=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_Xm=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_X0=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_Ym=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, max_iteration+1> gmres_V;
            Eigen::Matrix<double, max_iteration, max_iteration> gmres_Vm=Eigen::MatrixXd::Zero(max_iteration, max_iteration);
            Eigen::Matrix<double, max_iteration, 1> gmres_R0=Eigen::MatrixXd::Zero(max_iteration, 1);
            //gmres_X0を設定する
            gmres_X0=U;
            //初期残差gmres_R0を求める
            gmres_R0=calR0();
            gmres_V.col(0)=gmres_R0.normalized();
            Eigen::Matrix<double, max_iteration+1, 1> g;
            g(0, 0)=gmres_R0.norm();
            double h[max_iteration+1][max_iteration]{};
            Eigen::Matrix<double, max_iteration, max_iteration> Rm;
            double c[max_iteration]{};
            double s[max_iteration]{}; 
            Eigen::Matrix<double, u_size*N_step, 1> temp_sigma=Eigen::MatrixXd::Zero(u_size*N_step, 1);
            Eigen::Matrix<double, u_size*N_step, max_iteration+1> tempAv;
            Eigen::Matrix<double, u_size*N_step, 1> tempcalAv;
            for(int i=0; i<max_iteration; ++i){
                temp_sigma=Eigen::MatrixXd::Zero(u_size*N_step, 1);
                tempcalAv=calAv(gmres_V.col(i));
                for(int k=0; k<(i+1); ++k){
                    tempAv.col(k)=tempcalAv;
                    h[k][i]=tempAv.col(k).dot(gmres_V.col(k));
                    temp_sigma+=h[k][i]*gmres_V.col(k);
                }
                Eigen::Matrix<double, u_size*N_step, 1>temp_V=tempAv.col(i)-temp_sigma;
                //if(i==0)std::cout<<temp_V<<std::endl;
                h[i+1][i]=temp_V.norm();
                //if(i==0)std::cout<<h[1][0]<<std::endl;
                gmres_V.col(i+1)=temp_V/h[i+1][i];
                //if(i==0)std::cout<<gmres_V[1]<<std::endl;
                Rm(0, i)=h[0][i];
                for(int k=0; k<((i+1)-1); ++k){
                    double temp1=c[k]*Rm(k, i)+s[k]*h[k+1][i];
                    double temp2=-1*s[k]*Rm(k, i)+c[k]*h[k+1][i];
                    Rm(k, i)=temp1;
                    Rm(k+1, i)=temp2;
                }
                c[i]=Rm(i, i)/std::sqrt(Rm(i, i)*Rm(i, i)+h[i+1][i]*h[i+1][i]);
                s[i]=h[i+1][i]/std::sqrt(Rm(i, i)*Rm(i, i)+h[i+1][i]*h[i+1][i]);
                g(i+1, 0)=-1*s[i]*g(i, 0);
                g(i, 0)=c[i]*g(i, 0);
                Rm(i, i)=c[i]*Rm(i, i)+s[i]*h[i+1][i];
            }
            Eigen::Matrix<double, max_iteration, 1> Gm=g.block(0, 0, max_iteration, 1);
            gmres_Vm=gmres_V.block(0, 0, max_iteration, max_iteration);
            //後退代入によってRm*Ym=Gmを解く
            for(int i=(max_iteration-1); i>=0; --i){
                double temp_sigma_back=0;
                for(int k=(i+1); k<max_iteration; ++k){
                    temp_sigma_back+=Rm(i, k)*gmres_Ym[k];
                }
                gmres_Ym[i]=(Gm[i]-temp_sigma_back)/Rm(i, i);
            }
            gmres_Xm=gmres_X0+gmres_Vm*gmres_Ym;
            dU=gmres_Xm;
            U+=dU*ht;
            return U.block(0, 0, u_size, 1);
        }
        Eigen::Matrix<double, f_size*N_step, 1> calF(Eigen::Matrix<double, u_size*N_step, 1> _U, Eigen::Matrix<double, x_size, 1> _x){
            Eigen::Matrix<double, f_size*N_step, 1> F;
            //制約なし
            //0~Nまでx1, x2
            Eigen::Matrix<double, x_size*(N_step+1), 1> X_;
            //1~Nまでlamda1, lamda2
            Eigen::Matrix<double, x_size*N_step, 1> Lamda_;
            //x_(x*)を求める
            //x0*(t)=x(t)を計算する
            X_.block(0, 0, x_size, 1)=_x;
            //xi+1*=xi*+f(xi*,ui*,t+i*dtau)*dtauを計算する
            Eigen::Matrix<double, x_size, 1> prev_X_=_x;
            for(int i=1; i < (N_step+1); i++){
                X_.block(x_size*i, 0, x_size, 1)=prev_X_+calModel(prev_X_, _U.block(u_size*(i-1), 0, u_size, 1))*dt;
                prev_X_=X_.block(x_size*i, 0, x_size, 1);
            }
            //4
            //lamda_(lamda*)を求める
            //lamdaN*=(rphi/rx)^T(xN*,t+T)を計算する
            Lamda_.block(x_size*(N_step-1), 0, x_size, 1)=rphirx(X_.block(x_size*((N_step+1)-1), 0, x_size, 1)).transpose();
            //lamdai*=lamdai+1*+(rH/ru)^T*dtau
            Eigen::Matrix<double, x_size, 1> prev_Lamda_=Lamda_.block(x_size*(N_step-1), 0, x_size, 1);
            //逆順で解く
            //N_step-2の理由(N_step-1で最後のLamdaのグループなので(上でそこは計算してる),それの前だからN-2)
            for(int i=(N_step-2); i >= 0; --i){
                Lamda_.block(x_size*i, 0, x_size, 1)=prev_Lamda_+rHrx(X_.block(x_size*i, 0, x_size, 1), _U.block(u_size*i, 0, u_size, 1), prev_Lamda_).transpose()*dt;
                prev_Lamda_=Lamda_.block(x_size*i, 0, x_size, 1);
            }
            //Fを求める
            for(int i=0; i<N_step; i++){
                double lam_1=Lamda_((i*x_size), 0);
                double lam_2=Lamda_((i*x_size)+1, 0);
                double lam_3=Lamda_((i*x_size)+2, 0);
                double u_1=_U(i*u_size, 0);
                double u_2=_U((i*u_size+1), 0);
                double v_1=_U((i*u_size+2), 0);
                double v_2=_U((i*u_size+3), 0);
                double rho_1=_U((i*u_size+4), 0);
                double rho_2=_U((i*u_size+5), 0);
                double x_3=X_((i*x_size+2), 0);
                F((i*f_size),   0)=u_1+lam_1*std::cos(x_3)+lam_2*std::sin(x_3)+2*rho_1*u_1;
                F((i*f_size)+1, 0)=u_2+lam_3+2*rho_2*u_2;
                F((i*f_size)+2, 0)=-0.01+2*rho_1*v_1;
                F((i*f_size)+3, 0)=-0.01+2*rho_2*v_2;
                F((i*f_size)+4, 0)=u_1*u_1+v_1*v_1-1*1;
                F((i*f_size)+5, 0)=u_2*u_2+v_2*v_2-1.5*1.5;
            }
            return F;
        }   
        Eigen::Matrix<double, u_size*N_step, 1> calAv(Eigen::Matrix<double, u_size*N_step, 1> _V){
            Eigen::Matrix<double, u_size, 1> U0= U.block(0, 0, u_size,1);
            Eigen::Matrix<double, u_size*N_step, 1> Av=(calF(U+(_V*ht), X+calModel(X, U0)*ht)-calF(U, X+calModel(X, U0)*ht))/ht;
            return Av;
        }
        Eigen::Matrix<double, u_size*N_step, 1> calR0(){
            //U'(0)=U0を使用する
            Eigen::Matrix<double, x_size, 1> dX=calModel(X, U.block(0, 0, u_size,1))*ht;
            Eigen::Matrix<double, u_size*N_step, 1> R0=-1*zeta*calF(U, X)-(calF(U, X+dX)-calF(U, X))/ht-(calF(U+U*ht, X+dX)-calF(U, X+dX))/ht;
            return R0;
        }
        Eigen::Matrix<double, x_size, 1> X=Eigen::MatrixXd::Zero(x_size, 1);
    private:
        Eigen::Matrix<double, u_size*N_step, 1> U=Eigen::MatrixXd::Zero(u_size*N_step, 1);
        double dt;
        std::vector<double> save_x1;
        std::vector<double> save_x2;
        std::vector<double> save_x3;
        std::vector<double> save_u1;
        std::vector<double> save_u2;
};
int main(){
    constexpr double dt=0.01;
    constexpr int iteration_time=30;
    constexpr int iteration_num=iteration_time/dt;
    //現在の状態
    //x={x1, x2}
    Eigen::Matrix<double, x_size, 1> initX;
    //X(0)を測定する(初期値を代入する)
    initX << x1,
             x2,
             x3;
    NMPC nmpc(initX);
    std::vector<double> save_time;
    for(int i=1; i<iteration_num; ++i){
        double time=i*dt;
        save_time.push_back(time);
        Eigen::Matrix<double, u_size, 1> u=nmpc.CGMRES(time);
        nmpc.updateState(u, dt);
        std::cout<<i<<std::endl;
        //std::cout<<u<<std::endl;
    }
    while(1){
        nmpc.figGraph(save_time);
    }
}