#include<iostream>
#include<cmath>
#include<vector>
#include "/usr/include/eigen3/Eigen/Dense"
#include "/usr/include/eigen3/Eigen/Sparse"
#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/LU"

constexpr int x_size=2;
constexpr int u_size=3;
constexpr int f_size=3;
constexpr int raw_size=1;
constexpr int N_step=10;
constexpr double zeta=100.0;
constexpr double h=0.01;
constexpr double tf=1.0;
constexpr double alpha=0.5;
constexpr double error=0.001;

class SampleSystem{
    public:
        Eigen::Matrix<double, x_size, 1> x;
        SampleSystem(double _init_x1, double _init_x2){
            x(0, 0)=_init_x1;
            x(1, 0)=_init_x2;
        }
        void update_state(Eigen::Matrix<double, 2, 1> u){
            //ルンゲクッタを解く
        }
        double func_x_1(double y1, double y2, double u){
            double y_dot=y2;
            return y_dot;
        }
        double func_x_2(double y1, double y2, double u){
            double y_dot = (1. - y1*y1 - y2*y2) * y2 - y1 + u;
            return y_dot;
        }
};
class NMPCSimulatorSystem{
    public:
        NMPCSimulatorSystem(){
        }
        Eigen::Matrix<double, x_size*N_step, 1> cX_;
        Eigen::Matrix<double, x_size*N_step, 1> cLamda_;
        void calc_predict_and_adjoint_state(Eigen::Matrix<double, x_size, 1> x, Eigen::Matrix<double, u_size, 1> u, double dt){
            //制約なし
            //0~Nまでx1, x2, lamda1, lamda2
            Eigen::Matrix<double, x_size*N_step, 1> X_;
            Eigen::Matrix<double, x_size*N_step, 1> Lamda_;
            //x_(x*)を求める
            //x0*(t)=x(t)を計算する
            //部分行列を用いる
            X_.block(0, 0, x_size, 1)=_x;
            //xi+1*=xi*+f(xi*,ui*,t+i*dtau)*dtauを計算する
            Eigen::Matrix<double, x_size, 1> prev_X_=_x;
            for(int i=1; i < N_step; i++){
                X_.block(x_size*i, 0, x_size, 1)=prev_X_+calModel(prev_X_, _U.block(u_size*i, 0, u_size, 1), _t+i*dtau)*dtau;
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
                Eigen::Matrix<double, 1, x_size> temp_rHrx=rHrx(X_.block(x_size*i, 0, x_size, 1), _U.block(u_size*i, 0, u_size, 1), prev_Lamda_, _t+i*dtau);
                Lamda_.block(x_size*i, 0, x_size, 1)=prev_Lamda_+temp_rHrx.transpose()*dtau;
                prev_Lamda_=Lamda_.block(x_size*i, 0, x_size, 1);
            }
            cX_=X_;
            cLamda_=Lamda_;
        }
        Eigen::Matrix<double, x_size, 1> calModel(Eigen::Matrix<double, x_size, 1> x, double u){
            Eigen::Matrix<double, x_size, 1> x_dot;
            double x1=x(0, 0);
            double x2=x(1, 0);
            x_dot << x2,
                     (1-x1*x1-x2*x2)*x2-x1+u;
            return x_dot;
        }
};
class NMPCController_with_CGMRES{
    public:
        NMPCController_with_CGMRES(){
        }
        Eigen::Matrix<double, N_step, 1> calc_input(Eigen::Matrix<double, x_size, 1> x, double time){
            //calculating sampling time
            double dt = tf * (1.0 - std::exp(-alpha * time)) / float(N_step);

            //x_dot
            Eigen::Matrix<double, x_size, 1> x_dot=simulator.calModel(x, us[0]);

            Eigen::Matrix<double, x_size, 1> dx=x_dot*h;

            simulator.calc_predict_and_adjoint_state(x+dx, us, dt);
            
            Eigen::Matrix<double, x_size*N_step, 1> X_s=simulator.cX_;
            Eigen::Matrix<double, x_size*N_step, 1> Lamda_s=simulator.cLamda_;
        
            //Fxt
            Eigen::Matrix<double, f_size*N_step, 1> Fxt=_calc_f(X_s, Lamda_s, us, dummy_us, raws, dt);

            //F
            simulator.calc_predict_and_adjoint_state(x, us, dt);
            X_s=simulator.cX_;
            Lamda_s=simulator.cLamda_;

            Eigen::Matrix<double, f_size*N_step, 1> F=_calc_f(X_s, Lamda_s, us, dummy_us, raws, dt);

            Eigen::Matrix<double, f_size*N_step, 1> right = -zeta * F - ((Fxt - F) / h);

            Eigen::Matrix<double, N_step, 1> du = us*h;
            Eigen::Matrix<double, N_step, 1> ddummy_u=dummy_us*h;
            Eigen::Matrix<double, N_step, 1> draw=raws*h;

            simulator.calc_predict_and_adjoint_state(x+dx, us + du, dt);
            X_s=simulator.cX_;
            Lamda_s=simulator.cLamda_;

            Eigen::Matrix<double, f_size*N_step, 1> Fuxt=_calc_f(X_s, Lamda_s, us + du, self.dummy_us + ddummy_u,self.raws + draw, self.N, dt);

            Eigen::Matrix<double, f_size*N_step, 1> left=((Fuxt - Fxt) / self.ht);

            //calculationg cgmres
            Eigen::Matrix<double, f_size*N_step, 1> r0 = right - left;
            double r0_norm = r0.norm();
        
            vs = np.zeros((self.max_iteration, self.max_iteration + 1)) # 数×iterarion回数
        
            vs[:, 0] = r0 / r0_norm # 最初の基底を算出

            hs = np.zeros((self.max_iteration + 1, self.max_iteration + 1))

            e = np.zeros((self.max_iteration + 1, 1)) # in this case the state is 3(u and dummy_u)
            e[0] = 1.

            for i in range(self.max_iteration):
                du = vs[::3, i] * self.ht
                ddummy_u = vs[1::3, i] * self.ht
                draw = vs[2::3, i] * self.ht

            x_1s, x_2s, lam_1s, lam_2s = self.simulator.calc_predict_and_adjoint_state(x_1 + dx_1, x_2 + dx_2, self.us + du, self.N, dt)

            Fuxt = self._calc_f(x_1s, x_2s, lam_1s, lam_2s, self.us + du, self.dummy_us + ddummy_u,
                           self.raws + draw, self.N, dt)

            Av = (( Fuxt - Fxt) / self.ht)

            sum_Av = np.zeros(self.max_iteration)

            for j in range(i + 1): # グラムシュミットの直交化法です、和を取って差分を取って算出します
                hs[j, i] = np.dot(Av, vs[:, j])
                sum_Av = sum_Av + hs[j, i] * vs[:, j]

            v_est = Av - sum_Av

            hs[i+1, i] = np.linalg.norm(v_est)

            vs[:, i+1] = v_est / hs[i+1, i]

            inv_hs = np.linalg.pinv(hs[:i+1, :i]) # この辺は教科書（実時間の方）にのっています
            ys = np.dot(inv_hs, r0_norm * e[:i+1])

            judge_value = r0_norm * e[:i+1] - np.dot(hs[:i+1, :i], ys[:i])

            if np.linalg.norm(judge_value) < self.threshold or i == self.max_iteration-1:
                update_value = np.dot(vs[:, :i-1], ys_pre[:i-1]).flatten()
                du_new = du + update_value[::3]
                ddummy_u_new = ddummy_u + update_value[1::3]
                draw_new = draw + update_value[2::3]
                break

            ys_pre = ys
        
            //update
            self.us += du_new * self.ht
            self.dummy_us += ddummy_u_new * self.ht
            self.raws += draw_new * self.ht

            x_1s, x_2s, lam_1s, lam_2s = self.simulator.calc_predict_and_adjoint_state(x_1, x_2, self.us, self.N, dt)

            F = self._calc_f(x_1s, x_2s, lam_1s, lam_2s, self.us, self.dummy_us,
                            self.raws, self.N, dt)

            return us;
        }
        Eigen::Matrix<double, f_size*N_step, 1> _calc_f(Eigen::Matrix<double, x_size*N_step, 1> x, Eigen::Matrix<double, x_size*N_step, 1>lam, Eigen::Matrix<double, N_step, 1>us,  Eigen::Matrix<double, N_step, 1>dummy_us, Eigen::Matrix<double, raw_size*N_step, 1>raws, dt){
            Eigen::Matrix<double, f_size*N_step, 1> F;
            for(int i=0; i<N_step; i++){
                double lam_1=lam(i, 0);
                double lam_2=lam(i+1, 0);
                F((i*f_size),   0)=us(i, 0)+lam_2+2.0*raws(i, 0)*us(i, 0);
                F((i*f_size)+1, 0)=-0.01+2.0*raws(i, 0)*dummy_us(i, 0);
                F((i*f_size)+2, 0)=us(i, 0)*us(i, 0)+dummy_us(i, 0)*dummy_us(i, 0)-0.5*0.5;
            }
            return F;
        }
    private:
        NMPCSimulatorSystem simulator;

        Eigen::Matrix<double, N_step, 1> us;
        Eigen::Matrix<double, N_step, 1> dummy_us;
        Eigen::Matrix<double, N_step, 1> raws;
};

int main(){
    //simulation time
    double dt = 0.01;
    int iteration_time = 20.0;
    int iteration_num = int(iteration_time/dt);

    //plant
    double init_x_1=2.0;
    double init_x_2=0.0;
    SampleSystem plant_system(init_x_1, init_x_2);

    //controller
    NMPCController_with_CGMRES controller;

    for(int i=0; i<iteration_num; ++i){
        double time = float(i) * dt;
        Eigen::Matrix<double, x_size, 1> x=plant_system.x;
        Eigen::Matrix<double, N_step, 1>us=controller.calc_input(x, time);
        plant_system.update_state(us(0, 0));
    }
}