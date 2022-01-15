#include<iostream>
#include<cmath>
#include<vector>
#include"/usr/include/eigen3/Eigen/Dense"
#include"/usr/include/eigen3/Eigen/Sparse"
#include"/usr/include/eigen3/Eigen/Core"
#include"/usr/include/eigen3/Eigen/LU"
#include"/usr/include/eigen3/Eigen/QR"

int main(){
    Eigen::Matrix<double, 5, 5> a=Eigen::MatrixXd::Random(5, 5);
    std::cout<<a<<std::endl;
    std::cout<<a[:, :3]<<std::endl;
}