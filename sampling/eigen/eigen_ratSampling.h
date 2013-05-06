// This file is part of the KRYLSTAT function library
//
// Copyright (C) 2011 Erlend Aune <erlenda@math.ntnu.no>
//
// The KRYLSTAT library is free software; 
// you can redistribute it and/or modify it under the terms 
// of the GNU Lesser General Public License as published by
// the Free Software Foundation; either version 3 of the 
// License, or (at your option) any later version.
//
// Alternatively, you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of
// the License, or (at your option) any later version.
//
// The KRYLSTAT library is distributed in the 
// hope that it will be useful, but WITHOUT ANY WARRANTY; without
// even the implied warranty of MERCHANTABILITY or FITNESS FOR A
// PARTICULAR PURPOSE. See the GNU Lesser General Public License
// or the GNU General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License and a copy of the GNU General Public License along with
// the KRYLSTAT library. If not, see 
// <http://www.gnu.org/licenses/>.

#ifndef EIGEN_RATSAMPLING_H
#define EIGEN_RATSAMPLING_H

/*! Include openmp */
#include <omp.h>

/*! Including Eigen  */
#include <Eigen/Eigen>
//#include <unsupported/Eigen/SparseExtra>

/*! Including what I need from boost */
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

/*! Including the rest of the stuff */
//#include "../../cocg_m/basis/monitor.h"
#include "../../cocg_m/eigen/eigen_cocg_m.h"
#include "../../ratApps/sncndn.h"
#include "../../lanczos/eigen/eigen_lanczos.h"

using namespace Eigen;

typedef unsigned int int_type;


template <class LinearOperator, class Scalar, class Derived1>
void eigen_rat_sampling(const LinearOperator &A, const int_type nsamp, const MatrixBase<Derived1> &samples, 
			const Scalar &cgTol, const int_type &cgMaxIter);

template <class LinearOperator, class Scalar>
void eigen_ratSampling(const LinearOperator &A, const Scalar &min_eig, const Scalar &max_eig, 
		 const int_type nshifts, const int_type nsamp, 
		 Matrix<Scalar,Dynamic,Dynamic> &samples);







template <class LinearOperator, class Scalar, class Derived1>
void eigen_rat_sampling(const LinearOperator &A, const int_type nsamp, const MatrixBase<Derived1> &samples, 
			const Scalar &cgTol, const int_type &cgMaxIter)
{
  Scalar pi(3.141592653589793238462643383279L);
  Scalar minEig,maxEig;
  Matrix<Scalar,Dynamic,1> dummyX1(A.cols());
  Matrix<Scalar,Dynamic,1> dummyX2=Matrix<Scalar,Dynamic,1>::Zero(A.cols());
  
  boost::normal_distribution<Scalar> normDist(0,1);
  boost::mt11213b randGen;
  boost::variate_generator<boost::mt11213b, boost::normal_distribution<Scalar> > randNorm(randGen,normDist);

 #pragma omp parallel for 
  for (int_type jjj=0;jjj<A.rows();++jjj) {dummyX1(jjj)=randNorm();}

  default_monitor<Scalar> monEig(dummyX1, cgMaxIter, cgTol);
  
  
  
  eigen_lanc_extremal_eig_mult(A, dummyX1, dummyX2, minEig, maxEig, monEig);
  dummyX1=Matrix<Scalar,Dynamic,1>::Zero(1);
  dummyX2=Matrix<Scalar,Dynamic,1>::Zero(1);
  
  int nshifts(0);
  nshifts=static_cast<int>(-1.5*(std::log(maxEig/minEig)+3.)*std::log(cgTol)/(2.*pi*pi));
  std::cout << "nshifts=" << nshifts << "\n";
  
  Scalar* intConst = new Scalar;
  Scalar* wsq = new Scalar[nshifts];
  Scalar* dzdt = new Scalar[nshifts];
  
  sqrtIntPoints(nshifts,minEig,maxEig,intConst,wsq,dzdt);
  
  Matrix<Scalar,Dynamic,1> wsqEig(nshifts);
  
  // Quickfix. Wrapping should be better, but this is not performance intensive
  for (int_type iii=0;iii<nshifts;++iii){wsqEig(iii)=wsq[iii];}
  // std::cout << wsqEig << "\n";
  
  //samples.setZero();
  
#pragma omp parallel for
  for(int_type iii=0;iii<nsamp;++iii)
  {
    Matrix<Scalar,Dynamic,1> normalSample(A.cols());
    Matrix<Scalar,Dynamic,1> dummyX3(A.cols());
    Matrix<Scalar,Dynamic,Dynamic> x_sh(A.rows(),nshifts);
    for (int_type jjj=0;jjj<A.rows();++jjj)
    {
      normalSample(jjj)=randNorm();
    }
    default_monitor<Scalar> mon(normalSample);
    eigen_cocg_m(A,normalSample,dummyX3,x_sh,wsqEig,mon);
    for (int_type kkk=0;kkk<nshifts;++kkk)
    {
      const_cast< MatrixBase<Derived1>& >(samples).col(iii)=samples.col(iii)+dzdt[kkk]*x_sh.col(kkk);
    }
  }
  
  delete[] intConst;
  delete[] wsq;
  delete[] dzdt;
}







template <class LinearOperator, class Scalar>
void eigen_ratSampling(const LinearOperator &A, const Scalar &min_eig, const Scalar &max_eig, 
		 const int_type nshifts, const Scalar &cg_tol, const int_type nsamp, 
		 Matrix<Scalar,Dynamic,Dynamic> &samples)
{
  boost::normal_distribution<Scalar> normDist(0,1);
  boost::mt11213b randGen;
  boost::variate_generator<boost::mt11213b, boost::normal_distribution<Scalar> > randNorm(randGen,normDist);
  
  Scalar* intConst = new Scalar;
  Scalar* wsq = new Scalar[nshifts];
  Scalar* dzdt = new Scalar[nshifts];
  
  sqrtIntPoints(nshifts,min_eig,max_eig,intConst,wsq,dzdt);
  Matrix<Scalar,Dynamic,1> wsqEig(nshifts);
  
  // Quickfix. Wrapping should be better, but this is not performance intensive
  for (int_type iii=0;iii<nshifts;++iii){wsqEig(iii)=wsq[iii];}
  std::cout << wsqEig << "\n";
  
  // Ensure that samples are set to zero. Should be done beforehand...
  samples=0.*samples;
  
#pragma omp parallel for
  for(int_type iii=0;iii<nsamp;++iii)
  {
    Matrix<Scalar,Dynamic,1> normalSample(A.cols());
    Matrix<Scalar,Dynamic,1> dummyX(A.cols());
    Matrix<Scalar,Dynamic,Dynamic> x_sh(A.rows(),nshifts);
    for (int_type jjj=0;jjj<A.rows();++jjj)
    {
      normalSample(jjj)=randNorm();
    }
    default_monitor<Scalar> mon(normalSample);
    eigen_cocg_m(A,normalSample,dummyX,x_sh,wsqEig,mon);
    for (int_type kkk=0;kkk<nshifts;++kkk)
    {
      samples.col(iii)=samples.col(iii)+dzdt[kkk]*x_sh.col(kkk);
    }
  }
  
  
  delete[] intConst;
  delete[] wsq;
  delete[] dzdt;
  
}


#endif
