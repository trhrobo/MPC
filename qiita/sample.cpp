#include<iostream>
#include<cmath>
#include<vector>
#include "/usr/include/eigen3/Eigen/Dense"
#include "/usr/include/eigen3/Eigen/Sparse"
#include "/usr/include/eigen3/Eigen/Core"
#include "/usr/include/eigen3/Eigen/LU"

int main(){
    Eigen::Matrix<double, 3, 1> a;
    Eigen::Matrix<double, 3, 1> b;
    a<<1,2,3;
    b<<4,5,6;
    std::cout<<a.dot(b)<<std::endl;
}