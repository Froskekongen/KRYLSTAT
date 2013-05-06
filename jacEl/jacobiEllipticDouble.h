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

#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <stdio.h>
#include <iostream>
#include <complex>


/*! Include openmp */
#include <omp.h>
#include <typeinfo>

/*! Including Eigen  */
#include <Eigen/Eigen>
//#include <unsupported/Eigen/SparseExtra>

/*! Including limits */
#include <limits>

#ifndef JACOBI_ELLIPTIC_DOUBLE_H
#define JACOBI_ELLIPTIC_DOUBLE_H

using namespace Eigen;

namespace jacobiEllipticDouble
{
    template<class Scalar, class Derived, class IntegerType>
  void sqrtIntPoints(const IntegerType &N, const Scalar &minEig, const Scalar &maxEig, const Scalar &intConst, const MatrixBase<Derived> &wsq, const MatrixBase<Derived> &dzdt, int verbosity);
  
  template <class Scalar, class Derived, class IntegerType>
  void logIntPoints(const IntegerType &N, const Scalar &minEig, const Scalar &maxEig, Scalar &intConst,
		    const MatrixBase<Derived> &wsq, const MatrixBase<Derived> &dzdt, int verbosity);
  
  template<class Scalar>
  void ellKKP(const Scalar &L, Scalar &K, Scalar &Kp);
  
  template<class cpxDouble>
  cpxDouble poly_six(cpxDouble mmf);

  template<class cpxScalar,class Scalar>
  void recursiveSNCNDN(const cpxScalar &u, const cpxScalar &m, cpxScalar &sn, cpxScalar &cn, cpxScalar &dn, const Scalar dummy);

  template<class cpxScalar,class Scalar>
  void recursiveSNCNDN(const cpxScalar &u, const cpxScalar &m, cpxScalar &sn, cpxScalar &cn, cpxScalar &dn, const Scalar dummy)
  {
    const cpxScalar pi(3.141592653589793238462643383279L);
    const Scalar eps( std::numeric_limits<Scalar>::epsilon() );
    const Scalar CpLim=std::sqrt(eps);
    
    const cpxScalar zero(0.0L);
    const cpxScalar one(1.0L);
    const cpxScalar two(2.0L);
    const cpxScalar four(4.0L);
    const cpxScalar half(0.5L);
    const cpxScalar quart(0.25L);
    const cpxScalar three(3.0L);
    const cpxScalar five(5.0L);
    const cpxScalar six(6.0L);
    const cpxScalar cpOne(one);
    const Scalar mille(0.001L);
    
    cpxScalar mm=m;
    cpxScalar uu=u;
    cpxScalar mmf=mm/four;
    
    cpxScalar kappa,cpKapp,mu,v,sn1,cn1,dn1,denom,sinu,cosu;
    
    if( (mm.real()) < (four.real()*CpLim) )
    {
      sinu=sin(u);
      cosu=cos(u);
      sn=sinu+(mmf)*(sinu*cosu-u)*cosu;
      cn=cosu+(mmf)*(-sinu*cosu+u)*sinu;
      dn=one+(mmf)*(cosu*cosu-sinu*sinu-one);
    }
    else
    {
      
      if(mm.real()>mille)
      {
	kappa=(one-std::sqrt(one-mm))/(one+std::sqrt(one-mm));
      }
      else
      {
	kappa=poly_six(mmf);
      }
      cpKapp=cpxScalar(kappa);
    
      mu=std::pow(kappa,two);
      v=u/(cpOne+cpKapp);
    
      recursiveSNCNDN(v,mu,sn1,cn1,dn1,dummy);
    
      denom=cpOne+cpKapp*std::pow(sn1,two);
    
      sn=(cpOne+cpKapp)*sn1/denom;
      cn=cn1*dn1/denom;
      dn=(one-cpKapp*std::pow(sn1,two))/denom;
      
    }
    
  }

  template<class cpxDouble>
  cpxDouble poly_six(cpxDouble mmf)
  {
    const cpxDouble one(1.0L);
    const cpxDouble two(2.0L);
    const cpxDouble four(4.0L);
    const cpxDouble three(3.0L);
    const cpxDouble five(5.0L);
    const cpxDouble six(6.0L);
    cpxDouble kappa;
    kappa=cpxDouble(132.0L)*std::pow(mmf,six) + cpxDouble(42.0L)*std::pow(mmf,five) 
	+ cpxDouble(14.0L)*std::pow(mmf,four) + five*std::pow(mmf,three)
	+ two*std::pow(mmf,two) + one*std::pow(mmf,one);
    return kappa;
  }
  
  
  template<class Scalar, class Derived, class IntegerType>
  void logIntPoints(const IntegerType &N, const Scalar &minEig, const Scalar &maxEig, Scalar &intConst,
		    const MatrixBase<Derived> &wsq, const MatrixBase<Derived> &dzdt, int verbosity)
  {
    typedef typename Derived::Scalar cpxDouble;
    // typedef std::complex<Scalar> cpxDouble;
    const Scalar quart(0.25L);
    const Scalar one(1.0L);
    const Scalar eight(8.0L);
    const Scalar half(0.5L);
    const Scalar two(2.0L);
    const Scalar zero(0.0L);
    const Scalar pi(3.141592653589793238462643383279L);
    const Scalar realN(N);
    
    
    Scalar Mdiv=std::pow(maxEig/minEig,quart);
    Scalar Mmult=std::pow(maxEig*minEig,quart);
  
    Scalar k=(Mdiv-one)/(Mdiv+one);
    Scalar L=-std::log(k)/pi;
    
    Scalar K,Kp;
    ellKKP(L,K,Kp);
    
    intConst=-eight*K*Mmult/(k*pi*realN);
    
    //Derived t=Derived::Zero(N);
    Matrix<cpxDouble,Dynamic,1> t=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    for (int iii=0;iii<N;iii++)
    {
      cpxDouble temp1(zero,half*Kp),temp2(K,zero),temp3( (half+iii)*two*K/realN,zero );
      //t(iii)=cpxDouble(zero,half*Kp)-cpxDouble(K,zero) + cpxDouble((half+Scalar(iii))*two*K/realN,zero);
      t(iii)=temp1-temp2+temp3;
    }
    
    //Derived sn=Derived::Zero(N);
    //Derived cn=Derived::Zero(N);
    //Derived dn=Derived::Zero(N);
    Matrix<cpxDouble,Dynamic,1> sn=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    Matrix<cpxDouble,Dynamic,1> cn=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    Matrix<cpxDouble,Dynamic,1> dn=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    
    cpxDouble m=exp(-cpxDouble(two*pi*L));
    
    for (int iii=0;iii<N;++iii)
    {
      cpxDouble ttemp=t(iii),sntemp,cntemp,dntemp;
      Scalar dummy;
      recursiveSNCNDN(ttemp,m,sntemp,cntemp,dntemp,dummy);
      const_cast<cpxDouble &> (sn(iii))=sntemp;
      const_cast<cpxDouble &> (cn(iii))=cntemp;
      const_cast<cpxDouble &> (dn(iii))=dntemp;
      //sn(iii)=sntemp;
      //cn(iii)=cntemp;
      //dn(iii)=dntemp;
    }
    
    cpxDouble mmf=cpxDouble(Mmult);
    cpxDouble kres=cpxDouble(one)/cpxDouble(k);
    cpxDouble wTemp;
    cpxDouble wTempSQ;
    cpxDouble dzdtTemp1;
    cpxDouble dzdtTemp2;
    
    for (int iii=0;iii<N;++iii)
    {
      wTemp=mmf*(kres+sn(iii))/(kres-sn(iii));
      //wsq[iii]=mmf*(kres+sn[iii])/(kres-sn[iii]);
      wTempSQ=std::pow(wTemp,two);
      const_cast<cpxDouble &> (wsq(iii)) = cpxDouble(wTempSQ.real(),wTempSQ.imag());
    
      //wsq[iii]=std::pow(wsq[iii],2.);
    
      dzdtTemp1=cn(iii)*dn(iii)/std::pow(kres-sn(iii),two);
      dzdtTemp2=dzdtTemp1*std::log(wTempSQ)/wTemp;
    
      const_cast<cpxDouble &>(dzdt(iii)) = cpxDouble(dzdtTemp2.real(),dzdtTemp2.imag());
    
    //dzdt[iii]=cn[iii]*dn[iii]*std::pow(kres-sn[iii],2.);
    //dzdt[iii]=dzdt[iii]*std::log(wsq[iii])/std::sqrt(wsq[iii]);
    }
    const_cast<MatrixBase<Derived> &>(wsq)=-wsq;
    if(verbosity>2)
    {
      std::cout << "Shifts: " << wsq << "\n\n Contour: " << dzdt << "\n\n Integration constant: " << intConst << "\n\n\n";
    }
    
  }






  template<class Scalar>
  void ellKKP(const Scalar &L, Scalar &K, Scalar &Kp)
  {
    const Scalar pi(3.141592653589793238462643383279L);
    const Scalar eps( std::sqrt(std::numeric_limits<Scalar>::epsilon() ) );
    const Scalar two(2.0L);
    const Scalar one(1.0L);
    const Scalar zero(0.0L);
    const Scalar ten(10.0L);
    const Scalar four(4.0L);
  
    Scalar m=exp(-two*pi*(L));
    
    if(L>ten){K=pi/two;Kp=pi*L+std::log(four);return;}
  
    Scalar a0(one);
    Scalar b0=std::sqrt(one-m);
    Scalar s0=m;
    Scalar i1(zero);
    Scalar mm(one);
  
    Scalar a1(zero);
    Scalar b1(zero);
    Scalar c1(zero);
    Scalar w1(zero);
    
    
    while(std::abs(mm)>eps)
    {
      a1=(a0+b0)/two;
      b1=std::sqrt(a0*b0);
      c1=(a0-b0)/two;
      i1=i1+one;
      w1=std::pow(two,i1)*std::pow(c1,two);
      mm=w1;
      s0=s0+w1;
      a0=a1;
      b0=b1;
    }
    K=pi/(two*a1);
    
    if (  std::abs( one-m  )<eps ) // FIX THIS AT SOME POINT
    {
      K=(pi/(two*a1));
      // K->imag()=zero.real();
    }
    
    a0=one;
    b0=std::sqrt(m);
    s0=one-m;
    i1=zero;
    mm=one;
    
    while(std::abs(mm)>eps)
    {
      a1=(a0+b0)/two;
      b1=std::sqrt(a0*b0);
      c1=(a0-b0)/two;
      i1=i1+one;
      w1=std::pow(two,i1)*std::pow(c1,two);
      mm=w1;
      s0=s0+w1;
      a0=a1;
      b0=b1;
    }
    Kp=pi/(two*a1);
    
    if (  std::abs(m)<eps  )  // FIX THIS AT SOME POINT
    {
      Kp=(pi/(two*a1)); 
      //K->imag()=zero.real();
    }
    
    
  }
  
  
  template<class Scalar, class Derived, class IntegerType>
  void sqrtIntPoints(const IntegerType &N, const Scalar &minEig, const Scalar &maxEig, const Scalar &intConst, const MatrixBase<Derived> &wsq, const MatrixBase<Derived> &dzdt, int verbosity)
  {
    typedef std::complex<Scalar> cpxDouble;
  
    const Scalar pi(3.141592653589793238462643383279L);
  
    const Scalar k2(minEig/maxEig);
    const Scalar half(0.5L);
    const Scalar two(2.0L);
    const Scalar zero(0.0L);
    const Scalar realN(N);
    const Scalar L=-half*log(k2)/pi;
    
    Scalar K,Kp;
    ellKKP(L,K,Kp);

    const_cast<Scalar&>(intConst)=-two*Kp*sqrt(minEig)/(pi*realN);
    
    Matrix<cpxDouble,Dynamic,1> t(N);
    //mp_complex *t=new mp_complex[N];
    for(int iii=0;iii<N;iii++)
    {
      t(iii)=cpxDouble(zero,(half+static_cast<Scalar>(iii))*Kp/realN);
    }
    
    
//     mp_complex *sn=new mp_complex[N];
//     mp_complex *cn=new mp_complex[N];
//     mp_complex *dn=new mp_complex[N];
    Matrix<cpxDouble,Dynamic,1> sn=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    Matrix<cpxDouble,Dynamic,1> cn=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    Matrix<cpxDouble,Dynamic,1> dn=Matrix<cpxDouble,Dynamic,1>::Zero(N);
    
    cpxDouble m=exp(-two*pi*L);
    for(int iii=0;iii<N;iii++)
    {
      recursiveSNCNDN(t(iii),m,sn(iii),cn(iii),dn(iii),two);
    }
    
    
    //mp_complex minEigCpx(minEig);
    for(int iii=0;iii<N;iii++)
    {
      const_cast<Scalar&>(dzdt(iii))=(cn(iii)*dn(iii)).real();
      const_cast<Scalar&>(wsq(iii))=-(minEig*pow(sn(iii),two)).real();
    }
    
    if(verbosity>2)
    {
      std::cout << "intConst: " << intConst << "\n\n wsq: " << wsq << "\n\n dzdt: " << dzdt << "\n\n\n";
    }
  }
  
  
  
  
  
}




#endif