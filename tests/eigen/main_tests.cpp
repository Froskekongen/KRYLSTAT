#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET
#include <iostream>
#include <Eigen/Eigen>
#include <unsupported/Eigen/SparseExtra>

typedef unsigned int Integer;
typedef double ValueType;

using namespace Eigen;
int main()
{
  Integer sizeAmat(10000);
  DynamicSparseMatrix<ValueType,RowMajor,int> AB(sizeAmat,sizeAmat);
  AB.reserve(3*sizeAmat);
  for (Integer iii=0;iii<sizeAmat;++iii){AB.insert(iii,iii)=2.1;}
  for (Integer iii=0;iii<sizeAmat-1;++iii){AB.insert(iii+1,iii)=-1.;}
  for (Integer iii=0;iii<sizeAmat-1;++iii){AB.insert(iii,iii+1)=-1.;}
  AB.finalize();
  SparseMatrix<ValueType,RowMajor,Integer> A(AB);
  
  std::cout << "This works \n\n"
  
  
  SparseMatrix<ValueType,RowMajor,Integer> B=A*A;
  SparseMatrix<ValueType,RowMajor,int> C=A*A;
  
  std::cout << "This doesn't \n\n"
  
  
  
}