#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>
#include "eigen_ratSampling.h"
#include <boost/version.hpp>



int main()
{
  boost::mt11213b gen1;
  boost::normal_distribution<ValueType> ndnd(0,1);
  boost::variate_generator<boost::mt11213b, boost::normal_distribution<ValueType> > grog(gen1,ndnd);
  std::cout << "Boost version: " << BOOST_LIB_VERSION << "\n\n";
  
  int iii=0;
  while(iii<10)
  {
    std::cout << grog() << "\n";
    ++iii;
  }
  
  
  int_type sizeAmat=1000000;
  int_type nshifts=4;
  int_type Nsamples=10;
  
  DynamicSparseMatrix<double> AB(sizeAmat,sizeAmat);
  AB.reserve(3*sizeAmat);
  for (int iii=0;iii<sizeAmat;++iii){AB.insert(iii,iii)=2.001;}
  for (int iii=0;iii<sizeAmat-1;++iii){AB.insert(iii+1,iii)=-1.;}
  for (int iii=0;iii<sizeAmat-1;++iii){AB.insert(iii,iii+1)=-1.;}
  AB.finalize();
  SparseMatrix<double> A(AB);
  Matrix<double,Dynamic,Dynamic> samples(sizeAmat,Nsamples);
  
  double minEig=0.001,maxEig=4.5;
  double cg_tol=1e-3;
  int_type cgMaxIter=3000;
  
  eigen_ratSampling(A, minEig,maxEig,nshifts,cg_tol,Nsamples,samples);
  eigen_rat_sampling(A,Nsamples,samples,cg_tol,cgMaxIter);
  
  int *bog=A._innerIndexPtr();
  std::cout << "What is bog?  " << bog << "  What is its content?  " << *(bog+1) << "\n";
  
  return 1;
}