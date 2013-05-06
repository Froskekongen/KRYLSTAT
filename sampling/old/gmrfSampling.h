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

#ifndef GMRFSAMPLING_H
#define GMRFSAMPLING_H


#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>

#include <cusp/hyb_matrix.h>
#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/gallery/poisson.h>
#include <cusp/krylov/cg_m.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#include <jacEl/ellKKP.h>




typedef cusp::device_memory devMem;
typedef cusp::host_memory hostMem;

template<class VectorType> 
void arrayToCuspArray(int Nrat,double *intConst,double *wsq,double *dzdt,VectorType &cuspWsq,VectorType &cuspDzdt,VectorType &cuspIntConst)
{
  cuspIntConst[0]=*intConst;
  for (int iii=0;iii<Nrat;iii++)
  {
    cuspWsq[iii]=wsq[iii];
    cuspDzdt[iii]=dzdt[iii];
  } 
}



/*template <class LinearOperator,
          class VectorType1,
          class VectorType2,
          class VectorType3,
          class VectorType4,
          class Monitor
          class ValueType>
void cg_m_sampling(LinearOperator &A,ValueType &intConst,VectorType1 &wsq, VectorType1 &dzdt, VectorType2 &sampTemp,VectorType3 &outVec,VectorType4 &shiftedTemp. Monitor &mon)
{
  CUSP_PROFILE_SCOPED();
  sampTemp;
}*/



#endif
