#include<iostream>
#include<cmath>
#include<vector>
#include "/usr/include/eigen3/Eigen/Dense"
#include "/usr/include/eigen3/Eigen/Sparse"
#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/LU"

constexpr int x_size=2;
constexpr int N_step=10;
class SampleSystem{
    public:
        SampleSystem(double _init_x1, double _init_x2){
            x_1=_init_x1;
            x_2=_init_x2;
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
    private:
        double x_1;
        double x_2;
};
class NMPCSimulatorSystem{
    public:
        NMPCSimulatorSystem(){
        }
        void calc_predict_and_adjoint_state(){
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
            #if DEBUG0
            std::cout << "X_" << std::endl;
            std::cout << X_ << std::endl;
            #endif
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
        }
    private:
        Eigen::Matrix<double, 1, 1> x;
        Eigen::Matrix<double, 1, 1> lam;
        std::vector<double> x_1s;
        std::vector<double> x_2s;
        std::vector<double> lam_1s;
        std::vector<double> lam_2s;
};
class NMPCSimulatorSystem():
    def __init__(self):
        pass

    def calc_predict_and_adjoint_state(self, x_1, x_2, us, N, dt):
        x_1s, x_2s = self._calc_predict_states(x_1, x_2, us, N, dt) # by usin state equation
        lam_1s, lam_2s = self._calc_adjoint_states(x_1s, x_2s, us, N, dt) # by using adjoint equation

        return x_1s, x_2s, lam_1s, lam_2s

    def _calc_predict_states(self, x_1, x_2, us, N, dt):
        # initial state
        x_1s = [x_1]
        x_2s = [x_2]

        for i in range(N):
            temp_x_1, temp_x_2 = self._predict_state_with_oylar(x_1s[i], x_2s[i], us[i], dt)
            x_1s.append(temp_x_1)
            x_2s.append(temp_x_2)

        return x_1s, x_2s

    def _calc_adjoint_states(self, x_1s, x_2s, us, N, dt):
        # final state
        # final_state_func
        lam_1s = [x_1s[-1]]
        lam_2s = [x_2s[-1]]

        for i in range(N-1, 0, -1): 
            temp_lam_1, temp_lam_2 = self._adjoint_state_with_oylar(x_1s[i], x_2s[i], lam_1s[0] ,lam_2s[0], us[i], dt)
            lam_1s.insert(0, temp_lam_1)
            lam_2s.insert(0, temp_lam_2)

        return lam_1s, lam_2s

    def final_state_func(self):
        pass

    def _predict_state_with_oylar(self, x_1, x_2, u, dt):
        k0 = [0. for _ in range(2)]

        functions = [self.func_x_1, self.func_x_2]

        for i, func in enumerate(functions):
            k0[i] = dt * func(x_1, x_2, u)
                
        next_x_1 = x_1 + k0[0]
        next_x_2 = x_2 + k0[1]

        return next_x_1, next_x_2

    def func_x_1(self, y_1, y_2, u):
        y_dot = y_2
        return y_dot
    
    def func_x_2(self, y_1, y_2, u):
        y_dot = (1. - y_1**2 - y_2**2) * y_2 - y_1 + u
        return y_dot

    def _adjoint_state_with_oylar(self, x_1, x_2, lam_1, lam_2, u, dt):
        k0 = [0. for _ in range(2)]

        functions = [self._func_lam_1, self._func_lam_2]

        for i, func in enumerate(functions):
            k0[i] = dt * func(x_1, x_2, lam_1, lam_2, u)
                
        pre_lam_1 = lam_1 + k0[0]
        pre_lam_2 = lam_2 + k0[1]

        return pre_lam_1, pre_lam_2

    def _func_lam_1(self, y_1, y_2, y_3, y_4, u):
        y_dot = y_1 - (2. * y_1 * y_2 + 1.) * y_4
        return y_dot

    def _func_lam_2(self, y_1, y_2, y_3, y_4, u):
        y_dot = y_2 + y_3 + (-3. * (y_2**2) - y_1**2 + 1. ) * y_4
        return y_dot

class NMPCController_with_CGMRES():
    def __init__(self):
        """
        Parameters
        -----------
        None
        """
        # parameters
        self.zeta = 100. # 安定化ゲイン
        self.ht = 0.01 # 差分近似の幅
        self.tf = 1. # 最終時間
        self.alpha = 0.5 # 時間の上昇ゲイン
        self.N = 10 # 分割数
        self.threshold = 0.001 # break値

        self.input_num = 3 # dummy, 制約条件に対するuにも合わせた入力の数
        self.max_iteration = self.input_num * self.N

        # simulator
        self.simulator = NMPCSimulatorSystem()

        # initial
        self.us = np.zeros(self.N)
        self.dummy_us = np.ones(self.N) * 0.49
        self.raws = np.ones(self.N) * 0.011

        # for fig
        self.history_u = []
        self.history_dummy_u = []
        self.history_raw = []
        self.history_f = []

    def calc_input(self, x_1, x_2, time):
        # calculating sampling time
        dt = self.tf * (1. - np.exp(-self.alpha * time)) / float(self.N)

        # x_dot
        x_1_dot = self.simulator.func_x_1(x_1, x_2, self.us[0])
        x_2_dot = self.simulator.func_x_2(x_1, x_2, self.us[0])

        dx_1 = x_1_dot * self.ht
        dx_2 = x_2_dot * self.ht

        x_1s, x_2s, lam_1s, lam_2s = self.simulator.calc_predict_and_adjoint_state(x_1 + dx_1, x_2 + dx_2, self.us, self.N, dt)
        
        # Fxt
        Fxt = self._calc_f(x_1s, x_2s, lam_1s, lam_2s, self.us, self.dummy_us,
                            self.raws, self.N, dt)

        # F
        x_1s, x_2s, lam_1s, lam_2s = self.simulator.calc_predict_and_adjoint_state(x_1, x_2, self.us, self.N, dt)

        F = self._calc_f(x_1s, x_2s, lam_1s, lam_2s, self.us, self.dummy_us,
                            self.raws, self.N, dt)

        right = -self.zeta * F - ((Fxt - F) / self.ht)

        du = self.us * self.ht
        ddummy_u = self.dummy_us * self.ht
        draw = self.raws * self.ht

        x_1s, x_2s, lam_1s, lam_2s = self.simulator.calc_predict_and_adjoint_state(x_1 + dx_1, x_2 + dx_2, self.us + du, self.N, dt)

        Fuxt = self._calc_f(x_1s, x_2s, lam_1s, lam_2s, self.us + du, self.dummy_us + ddummy_u,
                           self.raws + draw, self.N, dt)

        left = ((Fuxt - Fxt) / self.ht)

        # calculationg cgmres
        r0 = right - left
        r0_norm = np.linalg.norm(r0)
        
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
        
        # update
        self.us += du_new * self.ht
        self.dummy_us += ddummy_u_new * self.ht
        self.raws += draw_new * self.ht

        x_1s, x_2s, lam_1s, lam_2s = self.simulator.calc_predict_and_adjoint_state(x_1, x_2, self.us, self.N, dt)

        F = self._calc_f(x_1s, x_2s, lam_1s, lam_2s, self.us, self.dummy_us,
                            self.raws, self.N, dt)

        print("check F = {0}".format(np.linalg.norm(F)))

        # for save
        self.history_f.append(np.linalg.norm(F))
        self.history_u.append(self.us[0])
        self.history_dummy_u.append(self.dummy_us[0])
        self.history_raw.append(self.raws[0])

        return self.us

    def _calc_f(self, x_1s, x_2s, lam_1s, lam_2s, us, dummy_us, raws, N, dt):
        F = []

        for i in range(N):
            F.append(us[i] + lam_2s[i] + 2. * raws[i] * us[i])
            F.append(-0.01 + 2. * raws[i] * dummy_us[i])
            F.append(us[i]**2 + dummy_us[i]**2 - 0.5**2)
        
        return np.array(F)

int main(){
    //simulation time
    double dt = 0.01
    int iteration_time = 20.
    int iteration_num = int(iteration_time/dt)

    //plant
    SampleSystem plant_system(init_x_1=2., init_x_2=0.)

    //controller
    controller = NMPCController_with_CGMRES()

    for(int i=0; i<iteration_num; ++i){
        double time = float(i) * dt
        double x_1 = plant_system.x_1
        double x_2 = plant_system.x_2
        # make inputコロナ
        us = controller.calc_input(x_1, x_2, time)
        # update state
        plant_system.update_state(us[0])
    }
}