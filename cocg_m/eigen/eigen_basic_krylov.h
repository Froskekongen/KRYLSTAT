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




#ifndef EIGEN_BASIC_KRYLOV_H
#define EIGEN_BASIC_KRYLOV_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>
#include <Eigen/Eigen>
#include <Eigen/Sparse>
//#include <unsupported/Eigen/SparseExtra>
#include "../basis/monitor.h"
#include <typeinfo>

using namespace Eigen;

namespace eigen_krylov{

template <class LinearOperator, class Monitor, class Derived>
void cg(const LinearOperator &A, const MatrixBase<Derived> &b, const MatrixBase<Derived> &x, Monitor &mon);

template <class LinearOperator, class LinearOperatorPrec, class Monitor, class Derived>
void cg(const LinearOperator &A, const LinearOperatorPrec &M, const MatrixBase<Derived> &b, const MatrixBase<Derived> &x, Monitor &mon);




template <class LinearOperator, class Monitor, class Derived>
void cg(const LinearOperator &A, const MatrixBase<Derived> &b, const MatrixBase<Derived> &x, Monitor &mon)
{
  typedef typename Derived::Scalar Scalar;
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());
  
  Derived r=b-A*x;
  Derived p=r;
  Derived ap(A.cols());
  
  Scalar alpha,beta;
  Scalar r_norm2=r.dot(r);
  Scalar r_norm2new;
  bool finish_flag=0;
  
  while(1)
  {
    ap=A*p;
    // r_norm2=mon.residual_norm()*mon.residual_norm();
    alpha=r_norm2/ap.dot(p);
    const_cast< MatrixBase<Derived>& >(x)=x+alpha*p;
    r=r-alpha*ap;
    
    if( mon.finished(r) )
    {
      break;
    } //Convergence check
    r_norm2new=mon.residual_norm()*mon.residual_norm();
    
    beta=r_norm2new/r_norm2;
    p=r+beta*p;
    r_norm2=r_norm2new;
    ++mon;
  }

}




template <class LinearOperator, class LinearOperatorPrec, class Monitor, class Derived>
void cg(const LinearOperator &A, const LinearOperatorPrec &M, const MatrixBase<Derived> &b, const MatrixBase<Derived> &x, Monitor &mon)
{
  typedef typename Derived::Scalar Scalar;
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());
  
  
  Derived r=b-A*x;
  Derived z=M*r;
  Derived p=z;
  Derived ap(A.cols());
  
  
  Scalar alpha,beta;
  Scalar rdotz=r.dot(z);
  
 
  
  while(!mon.finished(r))
  {
    ap=A*p;
    
    
    // std::cout << "laf";
    alpha=rdotz/ap.dot(p);
    
    const_cast< MatrixBase<Derived>& >(x)=x+alpha*p;
    r=r-alpha*ap;
    
    z=M*r;
    
    beta=r.dot(z)/rdotz;
    
    p=z+beta*p;
    
    rdotz=beta*rdotz;
    
    ++mon;
  }
}


template <class LinearOperator, class Monitor, class Derived>
void bicg_stab(const LinearOperator &A, const MatrixBase<Derived> &b, 
	       const MatrixBase<Derived> &x, Monitor &mon)
{
  typedef typename Derived::Scalar Scalar;
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());
  
  Derived r=b-A*x;
  Derived dualR=r;
  Derived p=r;
  Derived ap=Derived::Zero(A.cols());
  Derived s=Derived::Zero(A.cols());
  Derived as=Derived::Zero(A.cols());
  
  Scalar rDotDualR,rDotDnew,alpha,om,beta;
  
  while(!mon.finished(r))
  {
    ap=A*p;
    rDotDualR=r.dot(dualR);
    
    alpha=rDotDualR/ap.dot(dualR);
    s=r-alpha*ap;
    as=A*s;
    om=as.dot(s)/as.dot(as);
    const_cast< MatrixBase<Derived>& >(x)=x+alpha*p+om*s;
    r=s-om*as;
    
    rDotDnew=r.dot(dualR);
    beta=(rDotDnew/rDotDualR)*alpha/om;
    
    p=r+beta*(p-om*ap);
    ++mon;
  }
}



template <class LinearOperator, class LinearOperatorPrec, class Monitor, class Derived>
void bicg_stab(const LinearOperator &A, const LinearOperatorPrec &M, const MatrixBase<Derived> &b, 
	       const MatrixBase<Derived> &x, const MatrixBase<Derived> &dualR, Monitor &mon)
{
  typedef typename Derived::Scalar Scalar;
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows()); 
  
  Derived r=b-A*x;
  Derived p=r;
  Derived q(A.cols()), v(A.cols()), s(A.cols()), z(A.cols()), t(A.cols()), phi(A.cols());
  
  Scalar alpha, om, beta, rho_new;
  Scalar rho=dualR.dot(r);
  
  while(!mon.finished(r)) //something is weird with residual. Fix if can.
  {
    q=M*p;
    v=A*q;
    
    alpha=rho/dualR.dot(v);
    s=r-alpha*v;
    
    z=M*s;
    t=A*z;
    om=s.dot(t)/t.dot(t);
    
    r=s-om*t;
    const_cast< MatrixBase<Derived>& >(x)=x+alpha*q+om*z;
    
    rho_new=dualR.dot(r);
    beta=(alpha/om)*(rho_new/rho);
    
    p=r+beta*(p-om*v);
    
    rho=rho_new;
    ++mon;
  }

      
}






} // end namespace eigen_krylov

#endif
