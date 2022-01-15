#include<iostream>
#include<cmath>
#include<vector>
#include "/usr/include/eigen3/Eigen/Dense"
#include "/usr/include/eigen3/Eigen/Sparse"
#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/LU"
#include "matplotlibcpp.h"
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
            plt::named_plot("x1", _time, save_x1);
            plt::named_plot("x2", _time, save_x2);
            //plt::named_plot("u", _time, save_u);
            plt::legend();
            plt::pause(0.01);
        }
        void updateState(Eigen::Matrix<double, u_size, 1> _u){
            //状態Xを更新する
            double k0[2]{};
            double k1[2]{};
            double k2[2]{};
            double k3[2]{};

            double x1=X(0, 0);
            double x2=X(1, 0);
            double u=U(0, 0);
            save_x1.push_back(x1);
            save_x2.push_back(x2);
            save_u.push_back(u);
            for(int i=0; i<2; ++i){
                k0[i]=dt*Func(x1, x2, u, i);
            }
            for(int i=0; i<2; ++i){
                k1[i]=dt*Func(x1+k0[0]/2, x2+k0[1]/2, u, i);
            }
            for(int i=0; i<2; ++i){
                k2[i]=dt*Func(x1+k1[0]/2, x2+k1[1]/2, u, i);
            }
            for(int i=0; i<2; ++i){
                k3[i]=dt*Func(x1+k2[0]/2, x2+k2[1]/2, u, i);
            }
            x1+=(k0[0]+2*k1[0]+2*k2[0]+k3[0])/6;
            x2+=(k0[1]+2*k1[0]+2*k2[0]+k3[0])/6;
            X(0, 0)=x1;
            X(1, 0)=x2;
            /*X+=calModel(X, U.block(0, 0, u_size, 1))*dt;
            double x1=X(0, 0);
            double x2=X(1, 0);*/
            //std::cout<<x1<<" "<<x2<<std::endl;
        }
        double Func(double _x1, double _x2, double _u, double i){
            double ans{};
            if(i==0){
                ans=Func1(_x1, _x2, _u);
            }else{
                ans=Func2(_x1, _x2, _u);
            }
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
            dt=tf*(1-std::exp(-alpha*_time))/N_step;
            /*
            Eigen::Matrix<double, max_iteration, 1> gmres_R0=calR0();
            double r0_norm=gmres_R0.norm();
            Eigen::Matrix<double, max_iteration, max_iteration+1> Vs=Eigen::MatrixXd::Zero(max_iteration, max_iteration+1);
            Vs.col(0)=gmres_R0.normalized();
            Eigen::Matrix<double, max_iteration+1, max_iteration+1> Hs=Eigen::MatrixXd::Zero(max_iteration+1, max_iteration+1);
            Eigen::Matrix<double, max_iteration+1, 1> E;
            E(0, 0)=1;
            Eigen::Matrix<double, max_iteration, 1> dU_new;
            for(int i=0; i<max_iteration; ++i){
                Eigen::Matrix<double, max_iteration, 1> dU=Vs.col(i)*ht;
                Eigen::Matrix<double, u_size, 1> U0= U.block(0, 0, u_size,1);
                Eigen::Matrix<double, x_size, 1> dX=calModel(X, U0)*ht;
                Eigen::Matrix<double, max_iteration, 1> Av=((calF(X+dX, U+dU)-calF(X+dX, U))/ht);
                Eigen::Matrix<double, max_iteration, 1> sum_Av=Eigen::MatrixXd::Zero(max_iteration, 1);
                for(int j=0; j<(i+1); ++j){
                    Hs(j, i)=Av.dot(Vs.col(j));
                    sum_Av=sum_Av+Hs(j, i)*Vs.col(j);
                }
                Eigen::Matrix<double, max_iteration, 1> V_est=Av-sum_Av;
                Hs(i+1, i)=V_est.norm();
                Vs.col(i+1)=V_est/Hs(i+1, i);
                Eigen::Matrix<double, i+1, i> temp_Hs=Hs.block(0, 0, i+1, i);
                Eigen::Matrix<double, i+1, i> inv_hs=temp_Hs.completeOrthogonalDecomposition().pseudoInverse();
                Eigen::Matrix<double, i+1, 1> temp_E=E.block(i+1, 0);
                Eigen::Matrix<double, max_iteration, 1> ys=inv_hs.dot(r0_norm*temp_E);
                Eigen::Matrix<double, i+1, 1> judge_value=r0_norm*temp_E-temp_Hs*ys.block(i, 0);
                if(judge_value.norm()<threshold or i==(max_iteration-1)){
                    Eigen::Matrix<double, max_iteration, i-1> tempVs=Vs.dot(0, 0, max_iteration, i-1);
                    Eigen::Matrix<double, i-1, 1> temp_ys_pre=ys_pre.block(0, 0, i-1, 1);
                    Eigen::Matrix<double, max_iteration, 1> update_value=tempVs*temp_ys_pre;
                    dU_new+=update_value;
                    break;
                }
                Eigen::Matrix<double, max_iteration, 1> ys_pre=ys;
            }
            U+=dU_new*ht;*/
            //gmres法を用いてdUを求める
            Eigen::Matrix<double, max_iteration, 1> dU=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_Xm=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_X0=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_Ym=Eigen::MatrixXd::Zero(max_iteration, 1);
            Eigen::Matrix<double, max_iteration, 1> gmres_V[max_iteration];
            Eigen::Matrix<double, max_iteration, max_iteration> gmres_Vm=Eigen::MatrixXd::Zero(max_iteration, max_iteration);
            Eigen::Matrix<double, max_iteration, 1> gmres_R0=Eigen::MatrixXd::Zero(max_iteration, 1);
            //gmres_X0を設定する
            gmres_X0=U;
            //初期残差gmres_R0を求める
            gmres_R0=calR0();
            //std::cout<<calF(U, X)<<std::endl;
            //std::cout<<gmres_R0<<std::endl;
            gmres_V[0]=gmres_R0.normalized();
            double g[max_iteration+1]{};
            g[0]=gmres_R0.norm();
            double h[max_iteration+1][max_iteration]{};
            double r[max_iteration][max_iteration]{};
            double c[max_iteration]{};
            double s[max_iteration]{}; 
            for(int i=0; i<max_iteration; ++i){
                for(int k=0; k<(i+1); ++k){
                    Eigen::Matrix<double, u_size*N_step, 1> tempAv=calAv(gmres_V[i]);
                    //if(i==0)std::cout<<tempAv<<std::endl;
                    h[k][i]=tempAv.dot(gmres_V[k]);
                }
                Eigen::Matrix<double, u_size*N_step, 1> temp_sigma=Eigen::MatrixXd::Zero(u_size*N_step, 1);
                for(int k=0; k<(i+1); k++){
                    temp_sigma+=h[k][i]*gmres_V[k];
                }
                //if(i==0)std::cout<<temp_sigma<<std::endl;
                Eigen::Matrix<double, u_size*N_step, 1>temp_V=calAv(gmres_V[i])-temp_sigma;
                //if(i==0)std::cout<<temp_V<<std::endl;
                h[i+1][i]=temp_V.norm();
                //if(i==0)std::cout<<h[1][0]<<std::endl;
                gmres_V[i+1]=temp_V/h[i+1][i];
                //if(i==0)std::cout<<gmres_V[1]<<std::endl;
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
                //if(i==0)std::cout<<c[0]<<" "<<s[0]<<" "<<g[1]<<" "<<g[0]<<" "<<r[0]<<std::endl;
            }
            Eigen::Matrix<double, max_iteration, 1> Gm;
            for(int i=0; i<max_iteration; ++i){
                Gm(i, 0)=g[i];
            }
            Eigen::Matrix<double, max_iteration, max_iteration> Rm;
            for(int i=0; i<max_iteration; ++i){
                for(int k=0; k<max_iteration; ++k){
                    Rm(i, k)=r[i][k];
                }
            }
            //上で作ったgmres_V[i]をgmres_Vmに代入する
            for(int i=0; i<max_iteration; ++i){
                gmres_Vm.col(i)=gmres_V[i];
            }
            
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
            U+=dU*dt;
            return U.block(0, 0, u_size, 1);
        }
        int count=0;
        Eigen::Matrix<double, f_size*N_step, 1> calF(Eigen::Matrix<double, u_size*N_step, 1> _U, Eigen::Matrix<double, x_size, 1> _x){
            Eigen::Matrix<double, f_size*N_step, 1> F;
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
                X_.block(x_size*i, 0, x_size, 1)=prev_X_+calModel(prev_X_, _U.block(u_size*i, 0, u_size, 1))*dt;
                prev_X_=X_.block(x_size*i, 0, x_size, 1);
            }
            //4
            //lamda_(lamda*)を求める
            //lamdaN*=(rphi/rx)^T(xN*,t+T)を計算する
            Eigen::Matrix<double, 1, x_size> temp_rphirx=rphirx(X_.block(x_size*(N_step-1), 0, x_size, 1));
            Lamda_.block(x_size*(N_step-1), 0, x_size, 1)=temp_rphirx.transpose();
            //lamdai*=lamdai+1*+(rH/ru)^T*dtau
            Eigen::Matrix<double, x_size, 1> prev_Lamda_=Lamda_.block(x_size*(N_step-1), 0, x_size, 1);
            if(count==0){
                std::cout<<"prev_Lamda_"<<std::endl;
                std::cout<<prev_Lamda_<<std::endl;
            }
            //逆順で解く
            //N_step-2の理由(N_step-1で最後のLamdaのグループなので(上でそこは計算してる),それの前だからN-2)
            for(int i=(N_step-2); i >= 0; --i){
                Eigen::Matrix<double, 1, x_size> temp_rHrx=rHrx(X_.block(x_size*i, 0, x_size, 1), _U.block(u_size*i, 0, u_size, 1), prev_Lamda_);
                Lamda_.block(x_size*i, 0, x_size, 1)=prev_Lamda_+temp_rHrx.transpose()*dt;
                prev_Lamda_=Lamda_.block(x_size*i, 0, x_size, 1);
            }
            for(int i=0; i<N_step; i++){
                double lam_2=Lamda_((i*x_size)+1, 0);
                double u=_U(i*u_size, 0);
                double v=_U((i*u_size+1), 0);
                double rho=_U((i*u_size+2), 0);
                F((i*f_size),   0)=u+lam_2+2.0*rho*u;
                F((i*f_size)+1, 0)=-0.01+2.0*rho*v;
                F((i*f_size)+2, 0)=u*u+v*v-0.5*0.5;
            }
            if(count==0){
                std::cout<<"X_"<<std::endl;
                std::cout<<X_<<std::endl;
                std::cout<<"Lamda_"<<std::endl;
                std::cout<<Lamda_<<std::endl;
            }
            ++count;
            return F;
        }   
        Eigen::Matrix<double, u_size*N_step, 1> calAv(Eigen::Matrix<double, u_size*N_step, 1> _V){
            Eigen::Matrix<double, u_size, 1> U0= U.block(0, 0, u_size,1);
            Eigen::Matrix<double, x_size, 1> dX=calModel(X, U0)*ht;
            Eigen::Matrix<double, u_size*N_step, 1> temp1=calF(U+(_V*ht), X+dX);
            Eigen::Matrix<double, u_size*N_step, 1> temp2=calF(U, X+dX);
            Eigen::Matrix<double, u_size*N_step, 1> Av=(temp1-temp2)/ht;
            return Av;
        }
        Eigen::Matrix<double, u_size*N_step, 1> calR0(){
            //U'(0)=U0を使用する
            Eigen::Matrix<double, u_size*N_step, 1> dU=U*ht;
            Eigen::Matrix<double, u_size, 1> U0= U.block(0, 0, u_size,1);
            Eigen::Matrix<double, x_size, 1> dX=calModel(X, U0)*ht;
            Eigen::Matrix<double, u_size*N_step, 1> temp1=(calF(U, X+dX)-calF(U, X))/ht;
            Eigen::Matrix<double, u_size*N_step, 1> temp2=(calF(U+dU, X+dX)-calF(U, X+dX))/ht;
            Eigen::Matrix<double, u_size*N_step, 1> R0=-1*zeta*calF(U, X)-temp1-temp2;
            return R0;
        }
    private:
        Eigen::Matrix<double, x_size, 1> X=Eigen::MatrixXd::Zero(x_size, 1);
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
        //std::cout<<i<<std::endl;
        double time=i*dt;
        save_time.push_back(time);
        Eigen::Matrix<double, u_size, 1> u=nmpc.CGMRES(time);
        nmpc.updateState(u);
        //nmpc.figGraph(save_time);
    }
}