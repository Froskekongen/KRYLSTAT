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
#ifndef CG_H
#define CG_H



template <class LinearOperator,
          class Vector,
          class Monitor>
void cg(LinearOperator& A,
        Vector& x,
        Vector& b,Monitor &monitor);


template <class LinearOperator,
          class Vector,
          class Monitor,
          class Preconditioner>
void cg(LinearOperator& A,
        Vector& x,
        Vector& b,
        Monitor& monitor)
{
  assert(A.cols() == A.rows());
  assert(A.cols() == b.rows());  // sanity checks
  
  const int N=A.rows();
  
  Vector r=b;
  
  
  
}

#endif
