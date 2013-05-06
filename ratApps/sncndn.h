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

#ifndef SNCNDN_H
#define SNCNDN_H

#include <arprec/mp_real.h>
#include <arprec/mp_complex.h>

#include <math.h>
#include <complex>
#include <iostream>
#include <stdio.h>

typedef std::complex<double> cpxDbl;


void ellKKP(mp_real L,mp_real &K,mp_real &Kp);
void ellKKP(double L,double &K,double &Kp);

mp_complex poly_six(mp_complex mmf);
void recursiveSNCNDN(mp_complex u,mp_complex m, mp_complex &sn, mp_complex &cn, mp_complex &dn);
void recursiveSNCNDN(cpxDbl u,double m, cpxDbl &sn, cpxDbl &cn, cpxDbl &dn);
void recursiveSNCNDN(double u,double m, double &sn, double &cn, double &dn);

void sqrtIntPoints(int N,double minEig,double maxEig,double *intConst,double *wsq,double *dzdt);
void logIntPoints(int N, double minEig, double maxEig, double *intConst, cpxDbl *wsq, cpxDbl *dzdt);

#endif
