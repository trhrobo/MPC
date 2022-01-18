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
constexpr double x1=2.0;
constexpr double x2=0.0;

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
constexpr int max_iteration=u_size*N_step;
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
        void figGraph(std::vector<double> _time){
            plt::clf();
            //plt::plot(_time, save_x1);
            plt::plot(_time, save_x2);
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
            X+=calModel(X, _U.block(0, 0, u_size, 1))*dt;
            double x1=X(0, 0);
            double x2=X(1, 0);
            double u=_U(0, 0);
            save_x1.push_back(x1);
            save_x2.push_back(x2);
            save_u.push_back(u);
            //std::cout<<x1<<" "<<x2<<std::endl;
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
            //U={u, v, rho}
            double u=_U(0, 0);
            model_F << x2,
                       ((1-x1*x1-x2*x2)*x2-x1+u);
            return model_F;
        }
        Eigen::Matrix<double, 1, x_size> rphirx(Eigen::Matrix<double, x_size, 1> _X){
            Eigen::Matrix<double, 1, x_size> rphirx=_X.transpose();
            return rphirx;
        }
        Eigen::Matrix<double, 1, x_size> rHrx(Eigen::Matrix<double, x_size, 1> _x_,Eigen::Matrix<double, u_size, 1> _u, Eigen::Matrix<double, x_size, 1> _lamda_){
            double x1=_x_(0, 0);
            double x2=_x_(1, 0);
            double lamda1=_lamda_(0, 0);
            double lamda2=_lamda_(1, 0);
            Eigen::Matrix<double, 1, x_size> ans;
            ans<<x1-2*x1*x2*lamda2-lamda2, x2+lamda1+(-3*x2*x2-x1*x1+1)*lamda2;
            return ans;
        }
        Eigen::Matrix<double, u_size, 1> CGMRES(double _time){
            double time1=get_time_sec();
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
            double timep=get_time_sec();
            gmres_R0=calR0();
            double timeq=get_time_sec();
            //std::cout<<calF(U, X)<<std::endl;
            gmres_V.col(0)=gmres_R0.normalized();
            Eigen::Matrix<double, max_iteration+1, 1> g;
            g(0, 0)=gmres_R0.norm();
            double h[max_iteration+1][max_iteration]{};
            Eigen::Matrix<double, max_iteration, max_iteration> Rm;
            double c[max_iteration]{};
            double s[max_iteration]{}; 
            double time2=get_time_sec();
            Eigen::Matrix<double, u_size*N_step, 1> temp_sigma=Eigen::MatrixXd::Zero(u_size*N_step, 1);
            Eigen::Matrix<double, u_size*N_step, max_iteration+1> tempAv;
            Eigen::Matrix<double, u_size*N_step, 1> tempcalAv;
            for(int i=0; i<max_iteration; ++i){
                double timet=get_time_sec();
                temp_sigma=Eigen::MatrixXd::Zero(u_size*N_step, 1);
                tempcalAv=calAv(gmres_V.col(i));
                double timeab=get_time_sec();
                for(int k=0; k<(i+1); ++k){
                    double timeAv_start=get_time_sec();
                    tempAv.col(k)=tempcalAv;
                    double timeAv_end=get_time_sec();
                    //std::cout<<"Av"<<timeAv_end-timeAv_start<<std::endl;
                    //if(i==0)std::cout<<tempAv<<std::endl;
                    h[k][i]=tempAv.col(k).dot(gmres_V.col(k));
                    temp_sigma+=h[k][i]*gmres_V.col(k);
                }
                double timecd=get_time_sec();
                std::cout<<timecd-timeab<<std::endl;
                //if(i==0)std::cout<<temp_sigma<<std::endl;
                double timed=get_time_sec();
                //Eigen::Matrix<double, u_size*N_step, 1>tempAv=calAv(gmres_V.col(i));
                double timec=get_time_sec();
                Eigen::Matrix<double, u_size*N_step, 1>temp_V=tempAv.col(i)-temp_sigma;
                double timea=get_time_sec();
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
                //if(i==0)std::cout<<c[0]<<" "<<s[0]<<" "<<g[1]<<" "<<g[0]<<" "<<r[0]<<std::endl;
                double timeb=get_time_sec();
                //std::cout<<timeq-timep<<":"<<timeb-timed<<":"<<timea-timec<<":"<<timec-timed<<std::endl;
            }
            double time3=get_time_sec();
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
            double time4=get_time_sec();
            gmres_Xm=gmres_X0+gmres_Vm*gmres_Ym;
            dU=gmres_Xm;
            U+=dU*ht;
            //std::cout<<U.block(0, 0, u_size,1)<<std::endl;;
            std::cout<<time3-time2<<std::endl;
            return U.block(0, 0, u_size, 1);
        }
        Eigen::Matrix<double, f_size*N_step, 1> calF(Eigen::Matrix<double, u_size*N_step, 1> _U, Eigen::Matrix<double, x_size, 1> _x){
            double timef1=get_time_sec();
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
            double timef2=get_time_sec();
            Eigen::Matrix<double, x_size, 1> prev_X_=_x;
            for(int i=1; i < (N_step+1); i++){
                X_.block(x_size*i, 0, x_size, 1)=prev_X_+calModel(prev_X_, _U.block(u_size*(i-1), 0, u_size, 1))*dt;
                prev_X_=X_.block(x_size*i, 0, x_size, 1);
            }
            double timef3=get_time_sec();
            //4
            //lamda_(lamda*)を求める
            //lamdaN*=(rphi/rx)^T(xN*,t+T)を計算する
            Lamda_.block(x_size*(N_step-1), 0, x_size, 1)=rphirx(X_.block(x_size*((N_step+1)-1), 0, x_size, 1)).transpose();
            //lamdai*=lamdai+1*+(rH/ru)^T*dtau
            Eigen::Matrix<double, x_size, 1> prev_Lamda_=Lamda_.block(x_size*(N_step-1), 0, x_size, 1);
            //逆順で解く
            //N_step-2の理由(N_step-1で最後のLamdaのグループなので(上でそこは計算してる),それの前だからN-2)
            double timef4=get_time_sec();
            for(int i=(N_step-2); i >= 0; --i){
                Lamda_.block(x_size*i, 0, x_size, 1)=prev_Lamda_+rHrx(X_.block(x_size*i, 0, x_size, 1), _U.block(u_size*i, 0, u_size, 1), prev_Lamda_).transpose()*dt;
                prev_Lamda_=Lamda_.block(x_size*i, 0, x_size, 1);
            }
            double timef5=get_time_sec();
            //Fを求める
            for(int i=0; i<N_step; i++){
                double lam_2=Lamda_((i*x_size)+1, 0);
                double u=_U(i*u_size, 0);
                double v=_U((i*u_size+1), 0);
                double rho=_U((i*u_size+2), 0);
                F((i*f_size),   0)=u+lam_2+2.0*rho*u;
                F((i*f_size)+1, 0)=-0.01+2.0*rho*v;
                F((i*f_size)+2, 0)=u*u+v*v-0.5*0.5;
            }
            double timef6=get_time_sec();
            //std::cout<<timef2-timef1<<":"<<timef3-timef2<<":"<<timef4-timef3<<":"<<timef5-timef4<<timef6-timef5<<":"<<timef6-timef1<<std::endl;
            return F;
        }   
        Eigen::Matrix<double, u_size*N_step, 1> calAv(Eigen::Matrix<double, u_size*N_step, 1> _V){
            double timeAv1=get_time_sec();
            Eigen::Matrix<double, u_size, 1> U0= U.block(0, 0, u_size,1);
            Eigen::Matrix<double, u_size*N_step, 1> Av=(calF(U+(_V*ht), X+calModel(X, U0)*ht)-calF(U, X+calModel(X, U0)*ht))/ht;
            double timeAv2=get_time_sec();
            //std::cout<<timeAv2-timeAv1<<std::endl;
            return Av;
        }
        Eigen::Matrix<double, u_size*N_step, 1> calR0(){
            //U'(0)=U0を使用する
            double timecalr1=get_time_sec();
            Eigen::Matrix<double, x_size, 1> dX=calModel(X, U.block(0, 0, u_size,1))*ht;
            double timecalr2=get_time_sec();
            Eigen::Matrix<double, u_size*N_step, 1> R0=-1*zeta*calF(U, X)-(calF(U, X+dX)-calF(U, X))/ht-(calF(U+U*ht, X+dX)-calF(U, X+dX))/ht;
            double timecalr3=get_time_sec();
            //std::cout<<"calR0:"<<timecalr2-timecalr1<<":"<<timecalr3-timecalr2<<std::endl;
            return R0;
        }
        Eigen::Matrix<double, x_size, 1> X=Eigen::MatrixXd::Zero(x_size, 1);
        inline double get_time_sec(void){
            return static_cast<double>(duration_cast<nanoseconds>(steady_clock::now().time_since_epoch()).count())/1000000000;
        }
    private:
        //Eigen::Matrix<double, x_size, 1> X=Eigen::MatrixXd::Zero(x_size, 1);
        Eigen::Matrix<double, u_size*N_step, 1> U=Eigen::MatrixXd::Zero(u_size*N_step, 1);
        double dt;
        std::vector<double> save_x1;
        std::vector<double> save_x2;
        std::vector<double> save_u;
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
    std::vector<double> save_time;
    for(int i=1; i<iteration_num; ++i){
        double time=i*dt;
        save_time.push_back(time);
        double start=nmpc.get_time_sec();
        Eigen::Matrix<double, u_size, 1> u=nmpc.CGMRES(time);
        double end=nmpc.get_time_sec();
        //std::cout<<i<<":"<<end-start<<std::endl;
        nmpc.updateState(u, dt);
        nmpc.figGraph(save_time);
    }
}