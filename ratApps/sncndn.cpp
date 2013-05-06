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

#include "sncndn.h"

void ellKKP(mp_real L,mp_real &K,mp_real &Kp)
{
  const mp_real pi=mp_real::_pi;
  const mp_real eps=sqrt(mp_real::_eps);
  const mp_real two("2.0");
  const mp_real one("1.0");
  const mp_real zero("0.0");
  const mp_real ten("10.0");
  const mp_real four("4.0");
  
  mp_real m=exp(-two*pi*(L));
  
  if(L>ten){K=pi/two;Kp=pi*L+log(four);return;}
  
  mp_real a0(one);
  mp_real b0=sqrt(one-m);
  mp_real s0=m;
  mp_real i1(zero);
  mp_real mm(one);
  
  mp_real a1(zero);
  mp_real b1(zero);
  mp_real c1(zero);
  mp_real w1(zero);
  

  
  while(abs(mm)>eps)
  {
    a1=(a0+b0)/two;
    b1=sqrt(a0*b0);
    c1=(a0-b0)/two;
    i1=i1+one;
    w1=pow(two,i1)*pow(c1,two);
    mm=w1;
    s0=s0+w1;
    a0=a1;
    b0=b1;
  }
  K=pi/(two*a1);
  
  if (  abs( one-m  )<eps ) // FIX THIS AT SOME POINT
  {
    K=(pi/(two*a1));
    // K->imag()=zero.real();
  }
  
  a0=one;
  b0=sqrt(m);
  s0=one-m;
  i1=zero;
  mm=one;
  
  while(abs(mm)>eps)
  {
    a1=(a0+b0)/two;
    b1=sqrt(a0*b0);
    c1=(a0-b0)/two;
    i1=i1+one;
    w1=pow(two,i1)*pow(c1,two);
    mm=w1;
    s0=s0+w1;
    a0=a1;
    b0=b1;
  }
  Kp=pi/(two*a1);
  
  if (  abs(m)<eps  )  // FIX THIS AT SOME POINT
  {
    Kp=(pi/(two*a1)); 
    //K->imag()=zero.real();
  }
 
}




void ellKKP(double L,double &K,double &Kp)
{
  mp::mp_init(100,NULL,true);
  
  mp_real tempL=mp_real(L);
  mp_real tempK=mp_real(K);
  mp_real tempKp=mp_real(Kp);
  
  ellKKP(tempL,tempK,tempKp);
  
  K=dble(tempK);
  Kp=dble(tempKp);

  
  mp::mp_finalize();
}









mp_complex poly_six(mp_complex mmf)
{
  const mp_complex one("1.0");
  const mp_complex two("2.0");
  const mp_complex four("4.0");
  const mp_complex three("3.0");
  const mp_complex five("5.0");
  const mp_complex six("6.0");
  mp_complex kappa;
  kappa=mp_complex("132.0")*pow(mmf,six) + mp_complex("42.0")*pow(mmf,five) 
	+ mp_complex("14.0")*pow(mmf,four) + five*pow(mmf,three)
	+ two*pow(mmf,two) + one*pow(mmf,one);
  return kappa;
}

void recursiveSNCNDN(mp_complex u,mp_complex m, mp_complex &sn, mp_complex &cn, mp_complex &dn)
{
  const mp_complex pi=mp_complex(mp_real::_pi);
  const mp_real eps=mp_real::_eps;
  const mp_real CpLim=sqrt(eps);
  
  const mp_complex zero("0.0");
  const mp_complex one("1.0");
  const mp_complex two("2.0");
  const mp_complex four("4.0");
  const mp_complex half("0.5");
  const mp_complex quart("0.25");
  const mp_complex three("3.0");
  const mp_complex five("5.0");
  const mp_complex six("6.0");
  
  const mp_complex cpOne(one);
  
  const mp_real mille("0.001");
  
  mp_complex mm=m;
  mp_complex uu=u;
  
  mp_complex mmf=mm/four;
  
  
  mp_complex kappa;
  mp_complex cpKapp;
  mp_complex mu;
  mp_complex v;
  mp_complex sn1,cn1,dn1;
  mp_complex denom;
  
  mp_complex sinu;
  mp_complex cosu;
  
  if(mm.real<four.real*CpLim)
  {
    sinu=sin(u);
    cosu=cos(u);
    sn=sinu+(mmf)*(sinu*cosu-u)*cosu;
    cn=cosu+(mmf)*(-sinu*cosu+u)*sinu;
    dn=one+(mmf)*(cosu*cosu-sinu*sinu-one);
  }
  else
  {

    if (mm.real>mille)
    {
      kappa=(one-sqrt(one-mm))/(one+sqrt(one-mm));
    }
    else
    {
      kappa=poly_six(mmf);
    }
    cpKapp=mp_complex(kappa);
    
    mu=pow(kappa,two);
    v=u/(cpOne+cpKapp);
    
    recursiveSNCNDN(v,mu,sn1,cn1,dn1);
    
    denom=cpOne+cpKapp*pow(sn1,two);
    
    sn=(cpOne+cpKapp)*sn1/denom;
    cn=cn1*dn1/denom;
    dn=(one-cpKapp*pow(sn1,two))/denom;
  }
  
}


void recursiveSNCNDN(cpxDbl u,double m, cpxDbl &sn, cpxDbl &cn, cpxDbl &dn)
{
  mp::mp_init(100,NULL,true);
  
  mp_complex tempU=mp_complex(u.real(),u.imag());
  mp_complex tempM=mp_complex(m);
  mp_complex tempSN;
  mp_complex tempCN;
  mp_complex tempDN;
  
  recursiveSNCNDN(tempU,tempM,tempSN,tempCN,tempDN);
  
  sn=cpxDbl(dble(tempSN.real),dble(tempSN.imag));
  cn=cpxDbl(dble(tempCN.real),dble(tempCN.imag));
  dn=cpxDbl(dble(tempDN.real),dble(tempDN.imag));
  
  mp::mp_finalize();  
}


void recursiveSNCNDN(double u,double m, double &sn, double &cn, double &dn)
{
  mp::mp_init(100,NULL,true);
  
  mp_complex tempU=mp_complex(u);
  mp_complex tempM=mp_complex(m);
  mp_complex tempSN;
  mp_complex tempCN;
  mp_complex tempDN;
  
  recursiveSNCNDN(tempU,tempM,tempSN,tempCN,tempDN);
  
  sn=dble(tempSN.real);
  cn=dble(tempCN.real);
  dn=dble(tempDN.real);
  
  mp::mp_finalize();
}





void sqrtIntPoints(int N,double minEig,double maxEig,double *intConst,double *wsq,double *dzdt)
{
  mp::mp_init(100,NULL,true);
  
  const mp_real pi=mp_real::_pi;
  
  mp_real k2=mp_real(minEig)/mp_real(maxEig);
  mp_real L=-mp_real("0.5")*log(k2)/pi;
  mp_real K,Kp;
  ellKKP(L,K,Kp);
  
  mp_real half("0.5");
  mp_real two("2.0");
  mp_real zero("0.0");
  mp_real realN(N);
  
  mp_real tmpIntConst=-two*Kp*sqrt(mp_real(minEig))/(pi*realN);
  
  *intConst=dble(tmpIntConst);
  
  
  mp_complex *t=new mp_complex[N];
  for(int iii=0;iii<N;iii++)
  {
    
    t[iii]=mp_complex(zero,(half+mp_real(iii))*Kp/realN);
  }
  
  mp_complex *sn=new mp_complex[N];
  mp_complex *cn=new mp_complex[N];
  mp_complex *dn=new mp_complex[N];
  
  mp_complex m=exp(-mp_complex(two)*pi*L);
  for(int iii=0;iii<N;iii++)
  {
    recursiveSNCNDN(t[iii],m,sn[iii],cn[iii],dn[iii]);
  }
  
  mp_complex minEigCpx(minEig);
  for(int iii=0;iii<N;iii++)
  {
    dzdt[iii]=dble((cn[iii]*dn[iii]).real);
    wsq[iii]=dble(-(minEigCpx*pow(sn[iii],two)).real);
  }
  
  
  
  
  delete[] t;
  delete[] sn;
  delete[] cn;
  delete[] dn;
  
  mp::mp_finalize();
}



void logIntPoints(int N, double minEig, double maxEig, double *intConst, cpxDbl *wsq, cpxDbl *dzdt)
{
  mp::mp_init(100,NULL,true);
  
  
  const mp_real quart("0.25");
  const mp_real one("1.0");
  const mp_real eight("8.0");
  const mp_real half("0.5");
  const mp_real two("2.0");
  const mp_real zero("0.0");
  const mp_real pi=mp_real::_pi;
  const mp_real realN(N);
  
  mp_real Mdiv=pow(mp_real(maxEig)/mp_real(minEig),quart);
  mp_real Mmult=pow(mp_real(maxEig)*mp_real(minEig),quart);
  
  mp_real k=(Mdiv-one)/(Mdiv+one);
  mp_real L=-log(k)/pi;
  
  mp_real K,Kp;
  ellKKP(L,K,Kp);
  
  mp_real tmpIntConst=-eight*K*Mmult/(k*pi*realN);
  *intConst = dble(tmpIntConst);
  
  mp_complex *t=new mp_complex[N];
  for (int iii=0;iii<N;iii++)
  {
    t[iii]=mp_complex(zero,half*Kp)-mp_complex(K,zero) + mp_complex((half+mp_real(iii))*two*K/realN,zero);
  }
  
  mp_complex *sn=new mp_complex[N];
  mp_complex *cn=new mp_complex[N];
  mp_complex *dn=new mp_complex[N];
  
  
  mp_complex m=exp(-mp_complex(two*pi*L));
  
  for (int iii=0;iii<N;iii++)
  {
    recursiveSNCNDN(t[iii],m,sn[iii],cn[iii],dn[iii]);
  }
  
  
  mp_complex mmf=mp_complex(Mmult);
  mp_complex kres=mp_complex(one)/mp_complex(k);
  mp_complex wTemp;
  mp_complex wTempSQ;
  mp_complex dzdtTemp1;
  mp_complex dzdtTemp2;
  for (int iii=0;iii<N;iii++)
  {
    wTemp=mmf*(kres+sn[iii])/(kres-sn[iii]);
    //wsq[iii]=mmf*(kres+sn[iii])/(kres-sn[iii]);
    wTempSQ=pow(wTemp,two);
    wsq[iii]=cpxDbl(dble(wTempSQ.real),dble(wTempSQ.imag));
    
    //wsq[iii]=std::pow(wsq[iii],2.);
    
    dzdtTemp1=cn[iii]*dn[iii]/pow(kres-sn[iii],two);
    dzdtTemp2=dzdtTemp1*log(wTempSQ)/wTemp;
    
    dzdt[iii]=cpxDbl(dble(dzdtTemp2.real),dble(dzdtTemp2.imag));
    
    //dzdt[iii]=cn[iii]*dn[iii]*std::pow(kres-sn[iii],2.);
    //dzdt[iii]=dzdt[iii]*std::log(wsq[iii])/std::sqrt(wsq[iii]);
  }
  

  delete[] t;
  delete[] sn;
  delete[] cn;
  delete[] dn;
  mp::mp_finalize();
}
