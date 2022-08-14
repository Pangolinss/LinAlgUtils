#ifndef SYLVSOLVER_H
#define SYLVSOLVER_H

#include<eigen3/Eigen/Core>
#include<eigen3/Eigen/Eigenvalues>
#include<iostream>

//solve AX + XB = C
template<typename MatrixA, typename MatrixB, typename MatrixC>
class SylvSolverDense{
  using HessenA = Eigen::HessenbergDecomposition<MatrixA>;
  using HessenB = Eigen::HessenbergDecomposition<MatrixB>;
  //using VectorT = Eigen::Matrix<MatrixC::Scalar, Eigen::Dynamic, 1>;
  MatrixA& A;
  MatrixB& B;
  HessenA hessA;
  HessenB hessB;

  public:
  SylvSolverDense(MatrixA& _A, MatrixB& _B) : A{_A}, B{_B} {}
  Eigen::Matrix<typename MatrixA::Scalar, -1, -1> solve(const MatrixC& C)
  {
    Eigen::Matrix<typename MatrixA::Scalar , -1, -1> sol(A.rows(), B.cols());
    hessA.compute(A); //H = Q^TAQ
    hessB.compute(B.transpose()); //S = R^TBR 
    MatrixA H = hessA.matrixH();
    MatrixB S = hessB.matrixH();
    

    std::cout<<H<<std::endl<<std::endl;
    std::cout<<S<<std::endl<<std::endl;

    //transforms the equation into HY - YR^T = F
    MatrixC F = C;
    hessA.matrixQ().transpose().applyThisOnTheLeft(F); 
    hessB.matrixQ().applyThisOnTheRight(F);


    MatrixA HpskkI;// H + s_kk I
    Eigen::Matrix<typename MatrixA::Scalar, -1, -1> D(A.rows() * 2, A.cols() * 2); //matrix for the linear system we are solving in the else part
    int k = B.cols()-1;
    while(k >= 0)
    {
      if (k == 0 || S(k-1, k) == 0)
      {
        HpskkI = H;
        for (int i = 0; i < HpskkI.rows(); ++i)
          HpskkI(i,i) += S(k,k);

        auto fk = F.col(k);      

        for (int j = k+1; j < B.cols(); ++j)
        {
          fk -= S(k,j) * sol.col(j);
        }
        sol.col(k) = HpskkI.partialPivLu().solve(fk);

        k--;
      }
      else{
        D.setZero();
        D.block(0,0, A.rows(), A.cols()) = H;
        D.block(A.rows(), A.cols(), A.rows(), A.cols()) = H;

        D.block(0,0, A.rows(), A.cols()).diagonal().array() += S(k-1,k-1);
        D.block(0,A.cols(), A.rows(), A.cols()).diagonal().array() += S(k-1,k);
        D.block(A.rows(), 0, A.rows(), A.cols()).diagonal().array() += S(k,k-1);
        D.block(A.rows(), A.cols(), A.rows(), A.cols()).diagonal().array() += S(k,k);

        Eigen::Matrix<typename MatrixC::Scalar, -1, 2> fk1fk = F.block(0,k-1, F.rows(), 2);
        for (int j = k+1; j < B.cols(); ++j)
        {
          fk1fk.col(0) -= S(k-1,j) * sol.col(j);
          fk1fk.col(1) -= S(k,j) * sol.col(j);
        }
        Eigen::Map<Eigen::Matrix<typename MatrixA::Scalar, -1, 1>>(sol.block(0, k-1, sol.rows(), 2).data(), sol.rows()*2) = 
          D.partialPivLu().solve(Eigen::Map<Eigen::Matrix<typename MatrixC::Scalar, -1, 1>>(fk1fk.data(), fk1fk.size()));
        std::cout<<sol.block(0, k-1, sol.rows(), 2)<<std::endl<<std::endl;
        k-=2;
      }
    }

    hessA.matrixQ().applyThisOnTheLeft(sol);
    hessB.matrixQ().transpose().applyThisOnTheRight(sol);

    return sol;
  }
  

  public:

  
};

#endif