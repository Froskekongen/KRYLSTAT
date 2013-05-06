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

#ifndef COCG_M_H
#define COCG_M_H

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

//#include <Eigen/Eigen>
#include "monitor.h"
#include <typeinfo>
// using namespace Eigen;

typedef double ValueType;
typedef unsigned int int_type;
typedef std::complex<double> cpxDbl;

bool isnan(cpxDbl& a)
{
  return isnan(a.real()) || isnan(a.imag());
}


template <class LinearOperator,class Monitor, class ShiftMatrix, class ShiftVector, class Vector>
void cocg_m(LinearOperator &A, Vector &b, Vector &x, ShiftMatrix &x_sh, ShiftVector &shifts,Monitor &monitor);


template <class ShiftVector>
void first_shift(ShiftVector &z_old, ShiftVector &z_cur, ValueType &b_old, ValueType &b_cur, ValueType &alpha0,
		 ShiftVector &shifts, ShiftVector &z_new, ShiftVector &b_sh, int_type &len);

		 
template <class ShiftVector,class ShiftMatrix,class Vector>
void second_shift(ShiftMatrix &p_sh, ValueType &alpha0, ValueType &b_cur, ShiftVector &z_new, ShiftVector &z_cur,
		  ShiftVector &b_sh, Vector &r, int_type &len);


template <class LinearOperator,class Monitor, class ShiftMatrix, class ShiftVector, class Vector>
void cocg_m(LinearOperator &A, Vector &b, Vector &x, ShiftMatrix &x_sh, ShiftVector &shifts,Monitor &monitor)
{
  std::cout << "Initialising COCG-M \n";
  
  // Check complexity
  bool isComplex=0;
  if (typeid(shifts(0))==typeid(cpxDbl)){isComplex=1;};
  std::cout << "We have complex shifts: " << isComplex << "\n";
  
  
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());  // sanity checks
  
  //const int N=A.rows();
  const int_type ns=shifts.size(); //don't know if length exists. use size instead.
  int_type len = ns;
  
  // Nonshifted quantities
  Vector r=b;
  Vector p=r;
  // Vector x=0.*b; //in real version, x i excluded and only x_sh is evaluated
  ValueType b_norm=b.norm();
  
  ValueType b_old=1.;
  ValueType alpha0=0.;
  ValueType c_cur=b_norm*b_norm;
  
  // Shifted initialisers
  ShiftVector z_old=shifts;
  for (int_type iii=0;iii<ns;++iii){z_old(iii)=1;}
  //Matrix<ShiftScalar,Dynamic,1> z_old=MatrixXd::Ones(ns,1);
  ShiftVector z_cur=z_old;
  ShiftVector z_new=0.*z_old;
  ShiftVector b_sh=z_new;
  ShiftVector a_sh=z_new;
  
  x_sh=0*x_sh;
  ShiftMatrix p_sh=x_sh;

  if(!isComplex)
  {
    p_sh=b.replicate(1,ns);
  }
  else
  {
    //p_sh=b.replicate(1,ns).template cast<cpxDbl>();
    /*for (int_type iii=0; iii<A.cols(); ++iii)
    {
      p_sh.col(0)(iii)=cpxDbl(b(iii),0);
    }*/
    
    // Vector b1=b.template cast<cpxDbl>();
    /*for (int_type iii=0;iii<ns;++iii)
    {
      p_sh.col(iii)=b.template cast<Eigen::MatrixXcd>();
    }*/
  }
  
  // Noninitialised quantities
  Vector ap=b; // use 0.*b if strange
  ValueType b_cur;
  ValueType c_new;
  ShiftMatrix temp=x_sh;
  
  bool tempBool=0;
  
  
  std::cout << "CG is about to start \n";
  while (!monitor.finished(r))
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
    
    ++monitor;
  }

}


		 

template <class ShiftVector>
void first_shift(ShiftVector &z_old, ShiftVector &z_cur, ValueType &b_old, ValueType &b_cur, ValueType &alpha0,
		 ShiftVector &shifts, ShiftVector &z_new, ShiftVector &b_sh, int_type &len)
{
  //int_type ns=shifts.size();
  ShiftVector numer=0.*shifts;
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






#endif
