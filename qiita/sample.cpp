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
    Eigen::MatrixXd A(4, 5);
    //Eigen::Matrix<double, 3, 4> A;
    A(0,0)=1; A(0,1)=2;  A(0,2)=3; A(0,3)=4, A(0,4)=1;
    A(1,0)=5; A(1,1)=6;  A(1,2)=7; A(1,3)=8, A(1,4)=1;
    A(2,0)=9; A(2,1)=10;  A(2,2)=11; A(2,3)=12, A(2,4)=1;
    A(3,0)=9; A(3,1)=10;  A(3,2)=11; A(3,3)=12, A(3,4)=1;
    std::cout<<A<<std::endl;
    Eigen::MatrixXd pinv=A.block(0,0,3,4).completeOrthogonalDecomposition().pseudoInverse();
    //Eigen::Matrix<double, 4, 3> pinv=A.completeOrthogonalDecomposition().pseudoInverse();
    std::cout<<pinv<<std::endl;
}