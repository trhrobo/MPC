#include<iostream>
#include<cmath>
#include<vector>
#include"/usr/include/eigen3/Eigen/Dense"
#include"/usr/include/eigen3/Eigen/Sparse"
#include"/usr/include/eigen3/Eigen/Core"
#include"/usr/include/eigen3/Eigen/LU"
#include"/usr/include/eigen3/Eigen/QR"

using namespace std;
using namespace Eigen;

int main(){
    /*Eigen::Matrix<double, 3, 4> a;
    a<<1,2,3,4,
       5,6,7,8,
       9,10,11,12;
      // 10,11,12;
    Eigen::Matrix<double, 4, 3> ans;
    //ans=a.inverse();
    ans=a.completeOrthogonalDecomposition().pseudoInverse();
    std::cout<<a<<std::endl;
    std::cout<<ans<<std::endl;
    */
    Eigen::Matrix<double, 3, 4> A;
    A(0,0)=1; A(0,1)=2;  A(0,2)=3; A(0,3)=4;
    A(1,0)=5; A(1,1)=6;  A(1,2)=7; A(1,3)=8;
    A(2,0)=9; A(2,1)=10;  A(2,2)=11; A(2,3)=12;
    std::cout<<A<<std::endl;
    //Eigen::Matrix<double, 3, 4> pinv=A.completeOrthogonalDecomposition.pseudoInverse();
}