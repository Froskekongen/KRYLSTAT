#ifndef EIGEN_POW_H
#define EIGEN_POW_H


#include <Eigen/Eigen>
#include <Eigen/Sparse>
//#include <unsupported/Eigen/SparseExtra>
#include <iostream>

using namespace Eigen;

template <class Derived,class Int_Type>
void eigen_power(const MatrixBase<Derived> &A, Int_Type p, MatrixBase<Derived> &B);

template <class Derived1, class Derived2, class Int_Type>
void eigen_power_sparse(const SparseMatrixBase<Derived1> &A, Int_Type &p, SparseMatrixBase<Derived2> &B);





template <class Derived,class Int_Type>
void eigen_power(const MatrixBase<Derived> &A, Int_Type p, MatrixBase<Derived> &B)
{
  if (p > 0)
  {
    //Int_Type p = k;
    while (--p > 1)
    {
      B *= A;
    }
  }
  else
  {
    assert(false && "Not yet implemented");
    //return pow(Inverse(matrix), -power);
  }
}


template <class Derived1, class Derived2, class Int_Type>
void eigen_power_sparse(const SparseMatrixBase<Derived1> &A, const Int_Type &p,SparseMatrixBase<Derived2> &B)
{
  // B=A;
  if (p > 0)
  {
    
    for (Int_Type k = p;k>1;k=k-1)
    {
      B = A*B;
      // B.print();
      // std::cout << "Nonzeros: " << B.nonZeros() << "\n\n";
    }
  }
  else
  {
    assert(false && "Not yet implemented");
    //return pow(Inverse(matrix), -power);
  }
}




#endif
