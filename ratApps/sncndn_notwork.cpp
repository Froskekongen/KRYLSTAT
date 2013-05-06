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


void sncndn(double uu, double emmc, mp_real &sn, mp_real &cn, mp_real &dn) {
  static const mp_real CA=sqrt(mp_real::_eps); /* change to arbitrary precision */
  
  const mp_real zero("0.0");
  const mp_real one("1.0");
  const mp_real half("0.5");
  
  
  bool bo;
  int i,ii,l;
  
  const int numIter(13);
  
  mp_real a,b,c,d;
  
  mp_real emc;
  emc=mp_real(emmc);
  
  mp_real u;
  u=mp_real(uu);
  
  mp_real em[numIter],en[numIter];
  
  if (emc != zero) {
    bo=(emc < zero);
    if (bo) {
      d=one-emc;
      emc /= -one/d;
      u *= (d=sqrt(d));
    }
    a=one;
    dn=one;
    for (i=0;i<numIter;i++) {
      l=i;
      em[i]=a;
      en[i]=(emc=sqrt(emc));
      c=half*(a+emc);
      if (abs(a-emc) <= CA*a) break;
      emc *= a;
      a=c;
    }
    u *= c;
    sn=sin(u);
    cn=cos(u);
    if (sn != zero) {
      a=cn/sn;
      c *= a;
      for (ii=l;ii>=0;ii--) {
	b=em[ii];
	a *= c;
	c *= dn;
	dn=(en[ii]+a)/(b+a);
	a=c/b;
      }
      a=one/sqrt(c*c+one);
      sn=(sn >= zero ? a : -a);
      cn=c*sn;
    }
    if (bo) {
      a=dn;
      dn=cn;
      cn=a;
      sn /= d;
    }
  } else {
    cn=one/cosh(u);
    dn=cn;
    sn=tanh(u);
  }
}





void sncndn2(mp_real u,mp_real m,mp_real &sn,mp_real &cn,mp_real &dn,mp_real &am)
{
  // returns x = { sn, cn, dn, ph }
  const mp_real pi=mp_real::_pi;
  const mp_real eps=mp_real::_eps;
  const mp_real CpLim=sqrt(eps);
  const mp_real zero("0.0");
  const mp_real one("1.0");
  const mp_real two("2.0");
  const mp_real four("4.0");
  const mp_real half("0.5");
  const mp_real quart("0.25");
  

  mp_real ai, b, phi, t, twon;
  mp_real a[9];
  mp_real c[9];
  int i;

  // Check for special cases

  if ( m < zero || m > one )
  {
    // domain error
    // error( "Elliptic.Jacobi()", ERR_DOMAIN );
    printf( "Elliptic.Jacobi(). Domain error. m<0 or m>1");
    sn = zero;
    cn = zero;
    am = zero;
    dn = zero;
    return;
  }
  if( m < CpLim )
  {
    // m > 0 and near 0
    t = sin(u);
    b = cos(u);
    ai = quart * m * (u - t*b);
    sn = t - ai*b;
    cn = b + ai*t;
    am = u - ai;
    dn = one - half*m*t*t;
    return;
  }

  if( m >= one - CpLim )
  {
    // m < 1 and near 1
    ai = quart * (one-m);
    b = cosh(u);
    t = tanh(u);
    phi = one/b;
    twon = b * sinh(u);
    sn = t + ai * (twon - u)/(b*b);
    am = two*atan(exp(u)) - (half*pi) + ai*(twon - u)/b;
    ai *= t * phi;
    cn = phi - ai * (twon - u);
    dn = phi + ai * (twon + u);
    return;
  }

  // AGM scale
  a[0] = one;
  b = sqrt(one - m);
  c[0] = sqrt(m);
  twon = one;
  i = 0;

  while ( abs(c[i]/a[i]) > eps )
  {
    if( i > 7 )
    {
      //error( "Elliptic.Jacobi()", ERR_OVERFLOW );
      printf("Elliptic.Jacobi(), overflow error");

      
      // backward recurrence
      phi = twon * a[i] * u;
      do
      {
        t = c[i] * sin(phi) / a[i];
        b = phi;
        phi = half*(asin(t) + phi);
      }
      while ( --i != 0 );

      sn = sin(phi);
      t = cos(phi);
      cn = t;
      dn = t/cos(phi-b);
      am = phi;
      return;
    }
    ai = a[i];
    ++i;
    c[i] = half*( ai - b );
    t = sqrt( ai * b );
    a[i] = half*( ai + b );
    b = t;
    twon *= two;
  }

  // backward recurrence
  phi = twon * a[i] * u;
  do
  {
    t = c[i] * sin(phi) / a[i];
    b = phi;
    phi = half*(asin(t) + phi);
  }
  while ( --i != 0 );

  sn = sin(phi);
  t = cos(phi);
  cn = t;
  dn = t / cos(phi-b);
  am = phi;
  return;
}





void sncndn2(double u,double m,double &sn,double &cn,double &dn)
{
  mp::mp_init(100,NULL,true);
  
  mp_real tempU(u);
  mp_real tempM(m);
  mp_real tempSN(sn);
  mp_real tempCN(cn);
  mp_real tempDN(dn);
  mp_real tempAM;
  
  sncndn2(tempU,tempM,tempSN,tempCN,tempDN,tempAM);
  
  sn=dble(tempSN);
  cn=dble(tempCN);
  dn=dble(tempDN);
  
  mp::mp_finalize();
}



void sncndn2(mp_complex u,mp_real m, mp_complex &sn, mp_complex &cn, mp_complex &dn)
{
  mp_real x = u.real;
  mp_real y = u.imag;
  mp_real snx, cnx, dnx, amx;
  mp_real sn_, cn_, dn_, am_;
  mp_real one("1.0");

  // call real version
  sncndn2(x,m, snx, cnx, dnx, amx);
  sncndn2(y,one-m, sn_, cn_, dn_, am_);

  mp_real sny_i = sn_ / cn_;  // sny/i
  mp_real cny = one / cn_;
  mp_real dny = dn_ / cn_;

  mp_complex d = mp_complex(one) + mp_complex(m) * sqrt( snx * sny_i );
  sn = mp_complex( snx * cny * dny, sny_i * cnx * dnx ) / d;
  cn = mp_complex( cnx * cny, -snx * sny_i * dnx * dny ) / d;
  dn = mp_complex( dnx * dny, -m * snx * sny_i * cnx * cny ) / d;
}
