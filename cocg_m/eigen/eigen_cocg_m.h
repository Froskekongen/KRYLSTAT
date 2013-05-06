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


#ifndef EIGEN_COCG_M_H
#define EIGEN_COCG_M_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>
#include <Eigen/Eigen>
//#include <unsupported/Eigen/SparseExtra>
#include "../basis/monitor.h"
#include <typeinfo>

using namespace Eigen;

typedef double ValueType;
typedef unsigned int int_type;
typedef std::complex<double> cpxDbl;

typedef VectorXd Vd;
typedef MatrixXd Md;

typedef VectorXcd Vcpx;
typedef MatrixXcd Mcpx;



template <typename Scalar>
bool isnan(std::complex<Scalar> &a)
{
  return isnan(a.real()) || isnan(a.imag());
}




/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
template <class ShiftVector, class ShiftInput, typename Scal>
void first_shift(ShiftVector &z_old, ShiftVector &z_cur, Scal &b_old, Scal &b_cur, Scal &alpha0,
		 const ShiftInput &shifts, ShiftVector &z_new, ShiftVector &b_sh, int_type &len);

		 
template <class ShiftVector,class ShiftMatrix,class Vector>
void second_shift(ShiftMatrix &p_sh, ValueType &alpha0, ValueType &b_cur, ShiftVector &z_new, ShiftVector &z_cur,
		  ShiftVector &b_sh, Vector &r, int_type &len);
  
// template <class LinearOperator, class Monitor, class Scalar1, class Scalar2>
// void eigen_cocg_m(const LinearOperator &A, const Matrix<Scalar1,Dynamic,1> &b,
// 		   Matrix<Scalar1,Dynamic,1> &x, Matrix<Scalar2,Dynamic,Dynamic> &x_sh,
// 		   const Matrix<Scalar2,Dynamic,1> &shifts, Monitor &mon);

template <class LinearOperator, class Monitor, class Derived1, class Derived2, class Derived3>
void eigen_cocg_m(const LinearOperator &A, const MatrixBase<Derived1> &b,
		   const MatrixBase<Derived1> &x,const MatrixBase<Derived2> &x_sh,
		   const MatrixBase<Derived3> &shifts, Monitor &mon);

/*
template <class LinearOperator, class Monitor>
void eigen_cocg_m(const LinearOperator &A, const Vd &b, Vd &x, Mcpx &x_sh, const Vcpx &shifts, Monitor &mon);


// double overloaded
template <class LinearOperator, class Monitor>
void eigen_cocg_m(const LinearOperator &A, const Vd &b, Vd &x, Md &x_sh, const Vd &shifts, Monitor &mon);
*/
/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/



/*template <class LinearOperator, class Monitor, class Scalar1, class Scalar2>
void eigen_cocg_m(const LinearOperator &A, const Matrix<Scalar1,Dynamic,1> &b,
		   Matrix<Scalar1,Dynamic,1> &x, Matrix<Scalar2,Dynamic,Dynamic> &x_sh,
		   const Matrix<Scalar2,Dynamic,1> &shifts, Monitor &mon)
{
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());
  
  int_type N=A.cols();
  
  const int_type ns=shifts.size();
  int_type len=ns;
  
  // Nonshifted quantities
  Matrix<Scalar1,Dynamic,1> r=b;
  Matrix<Scalar1,Dynamic,1> p=r;
  
  Scalar1 b_norm=b.norm();
  Scalar1 b_old=1.;
  Scalar1 alpha0=0;
  Scalar1 c_cur=b_norm*b_norm;
  
  // Shifted initialisers
  Matrix<Scalar2,Dynamic,1> z_old=Matrix<Scalar2,Dynamic,1>::Ones(ns);
  Matrix<Scalar2,Dynamic,1> z_cur=Matrix<Scalar2,Dynamic,1>::Ones(ns);
  Matrix<Scalar2,Dynamic,1> z_new=Matrix<Scalar2,Dynamic,1>::Zero(ns);
  Matrix<Scalar2,Dynamic,1> b_sh=Matrix<Scalar2,Dynamic,1>::Zero(ns);
  // Vcpx a_sh=Vcpx::Zero(ns);
  
  x_sh=0.*x_sh;
  
  // This is the trouble maker
  Matrix<Scalar2,Dynamic,Dynamic> p_sh=b.replicate(1,ns).template cast<Scalar2>();

  
  // Noninitialised quantities
  Matrix<Scalar1,Dynamic,1> ap(N);
  Scalar1 b_cur;
  Scalar1 c_new;
  Matrix<Scalar2,Dynamic,Dynamic> temp(N,ns);
  
  bool tempBool=0;
  
  std::cout << "CG is about to start \n";
  while (!mon.finished(r))
  {
    ap=A*p;
    b_cur=-c_cur/(p.transpose()*ap);
    x=x-b_cur*p;
    
    // Compute shifted stuff
    first_shift(z_old,z_cur,b_old,b_cur,alpha0,shifts,z_new,b_sh,len);
    
    for(int_type jjj=0;jjj<len;++jjj)
    {
      temp.col(jjj)=b_sh(jjj)*p_sh.col(jjj);
      for(int_type kkk=0;kkk<b.size();++kkk)
      {
	if (isnan(temp.col(jjj)(kkk)))
	{
	  len=len-1;
	  tempBool=1;
	  break;
	}
      }
      if (tempBool){break;}
      x_sh.col(jjj)=x_sh.col(jjj)-b_sh(jjj)*p_sh.col(jjj);
    }
    // end compute shifted stuff
    
    r=r+b_cur*ap;
    c_new=r.transpose()*r;
    alpha0=c_new/c_cur;
    p=r+alpha0*p;
    
    // Compute shifted stuff
    second_shift(p_sh, alpha0, b_cur, z_new, z_cur, b_sh, r, len);
    // end compute shifted stuff
    
    z_old=z_cur; z_cur=z_new;
    b_old=b_cur; c_cur=c_new;
    
    ++mon;
  }
}*/


template <class LinearOperator, class Monitor, class Derived1, class Derived2, class Derived3>
void eigen_cocg_m(const LinearOperator &A, const MatrixBase<Derived1> &b,
		   const MatrixBase<Derived1> &x,const MatrixBase<Derived2> &x_sh,
		   const MatrixBase<Derived3> &shifts, Monitor &mon)
{
  typedef typename Derived1::Scalar Scalar1;
  typedef typename Derived2::Scalar Scalar2;
  //typedef typename internal::plain_column_type<Derived1> ColType1;
  
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());
  
  int_type N=A.cols();
  
  const int_type ns=shifts.size();
  int_type len=ns;
  
  // Nonshifted quantities
  Derived1 r=b;
  Derived1 p=r;
  
  Scalar1 b_norm=b.norm();
  Scalar1 b_old=1.;
  Scalar1 alpha0=0;
  Scalar1 c_cur=b_norm*b_norm;
  
  
  // Shifted initialisers
  Derived3 z_old(Derived3::Ones(ns));
  Derived3 z_cur(Derived3::Ones(ns));
  Derived3 z_new(Derived3::Zero(ns));
  Derived3 b_sh(Derived3::Zero(ns));
  // Vcpx a_sh=Vcpx::Zero(ns);
  
  //x_sh=0.*x_sh;
  // const_cast< MatrixBase<Derived2>& >(x_sh)=Scalar2(0.)*x_sh;
  
  // This is the trouble maker
  Derived2 p_sh=b.replicate(1,ns).template cast<Scalar2>();

  
  // Noninitialised quantities
  Derived1 ap(N);
  Scalar1 b_cur;
  Scalar1 c_new;
  Derived2 temp(N,ns);
  
  bool tempBool=0;
  
  // std::cout << "CG is about to start \n";
  while (!mon.finished(r))
  {
    ap=A*p;
    b_cur=-c_cur/(p.transpose()*ap);
    const_cast< MatrixBase<Derived1>& >(x)=x-b_cur*p;
    
    // Compute shifted stuff
    first_shift(z_old,z_cur,b_old,b_cur,alpha0,shifts,z_new,b_sh,len);
    
    for(int_type jjj=0;jjj<len;++jjj)
    {
      temp.col(jjj)=b_sh(jjj)*p_sh.col(jjj);
      /*for(int_type kkk=0;kkk<b.size();++kkk)
      {
	if (isnan(temp.col(jjj)(kkk)))
	{
	  len=len-1;
	  tempBool=1;
	  break;
	}
      }*/
      if ( (temp.col(jjj).array()!=temp.col(jjj).array()).any() )
      {
	len=len-1;
	tempBool=1;
      }
      
      if (tempBool){break;}
      const_cast< MatrixBase<Derived2>& >(x_sh).col(jjj)=x_sh.col(jjj)-b_sh(jjj)*p_sh.col(jjj);
    }
    // end compute shifted stuff
    
    r=r+b_cur*ap;
    c_new=r.transpose()*r;
    alpha0=c_new/c_cur;
    p=r+alpha0*p;
    
    // Compute shifted stuff
    second_shift(p_sh, alpha0, b_cur, z_new, z_cur, b_sh, r, len);
    // end compute shifted stuff
    
    z_old=z_cur; z_cur=z_new;
    b_old=b_cur; c_cur=c_new;
    tempBool=0;
    
    ++mon;
  }
}




/*! Shifted parameter templates */
template <class ShiftVector, class ShiftInput, typename Scal>
void first_shift(ShiftVector &z_old, ShiftVector &z_cur, Scal &b_old, Scal &b_cur, Scal &alpha0,
		 const ShiftInput &shifts, ShiftVector &z_new, ShiftVector &b_sh, int_type &len)
{
  //int_type ns=shifts.size();
  ShiftVector numer=ShiftVector::Zero(len);
  ShiftVector denom=numer;

  for (int_type iii=0;iii<len;++iii)
  {
    numer(iii)=z_old(iii)*z_cur(iii)*b_old;
    denom(iii)=b_cur*alpha0*(z_old(iii)-z_cur(iii));
    denom(iii)=denom(iii)+b_old*z_old(iii)*(1.-b_cur*shifts(iii));
    z_new(iii)=numer(iii)/denom(iii);
    b_sh(iii)=b_cur*z_new(iii)/z_cur(iii);

    if (isnan(b_sh(iii)))
    {
      len=len-1;
      return;
    }    
    
  }
}


template <class ShiftVector,class ShiftMatrix,class Vector>
void second_shift(ShiftMatrix &p_sh, ValueType &alpha0, ValueType &b_cur, ShiftVector &z_new, ShiftVector &z_cur,
		  ShiftVector &b_sh, Vector &r, int_type &len)
{
  ShiftVector a_sh(len);
  ValueType c1=alpha0/b_cur;
  ShiftVector numer(len);
  for (int_type iii=0;iii<len;++iii)
  {
    numer(iii)=c1*z_new(iii)*b_sh(iii);
    a_sh(iii)=numer(iii)/z_cur(iii);
    if (isnan(a_sh(iii)))
    {
      len=len-1;
      return;
    }
    p_sh.col(iii)=z_new(iii)*r+a_sh(iii)*p_sh.col(iii);
  }
}
/*! Shifted parameter templates */



#endif
