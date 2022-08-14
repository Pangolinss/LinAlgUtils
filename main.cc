#include<iostream>
#include<eigen3/Eigen/Core>
#include"SylvSolver.h"

template<typename Scalar>
using MatrixT = Eigen::Matrix<Scalar, -1, -1>;

int main()
{
  Eigen::Matrix<double, -1, -1> A(5,5);
  Eigen::Matrix<double, -1, -1> B(5,5);
  Eigen::Matrix<double, -1, -1> C(5,5);
  B.setIdentity();
  A.setZero();
  C.setZero();
  for (int i = 0; i < 5; i++)
  {
    A(4-i,i) = 2;
  }
  A(2,2) = 3;
  for (int i = 0; i < 5; i++)
    for (int j = 0; j < 5; j++)
      C(i,j) = 2.0f;
  SylvSolverDense<MatrixT<double>, MatrixT<double>, MatrixT<double>> sylv(B, A);
  MatrixT<double> sol = sylv.solve(C);
  std::cerr<<B * sol + sol * A - C;
  return 0;
}