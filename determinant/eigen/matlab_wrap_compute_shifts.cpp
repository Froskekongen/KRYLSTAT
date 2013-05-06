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

#include <iostream>
#include "../../jacEl/jacobiEllipticDouble.h"
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
    {
        mexErrMsgTxt("Need 4 input arguments: minEig,maxEig,nshifts,verbosity");
    }
    if (nlhs>3)
    {
        mexErrMsgTxt("Maximum of 3 output arguments");
    }
    
    double minEig(mxGetScalar(prhs[0])),maxEig(mxGetScalar(prhs[1]));
    int nshifts(mxGetScalar(prhs[2])),verbosity(mxGetScalar(prhs[3]));
    
    plhs[0]=mxCreateDoubleMatrix(nshifts,1,mxCOMPLEX);
    plhs[1]=mxCreateDoubleMatrix(nshifts,1,mxCOMPLEX);
    plhs[2]=mxCreateDoubleMatrix(1,1,mxREAL);
    double *wsqPtrReal=mxGetPr(plhs[0]);
    double *wsqPtrImag=mxGetPi(plhs[0]);
    double *dzdtPtrReal=mxGetPr(plhs[1]);
    double *dzdtPtrImag=mxGetPi(plhs[1]);
    
    double *intConstPtr;
    intConstPtr=mxGetPr(plhs[2]);
    double intConst;

    VectorXcd wsq=VectorXcd::Zero(nshifts),dzdt=VectorXcd::Zero(nshifts);
    //Map<VectorXcd> wsqEig(wsqPtr,nshifts),dzdtEig(dzdtPtr,nshifts);
    //VectorXcd a(wsqEig),b(dzdtEig);
    jacobiEllipticDouble::logIntPoints(nshifts, minEig, maxEig, intConst, wsq, dzdt, verbosity);
    *intConstPtr=intConst;
    for (std::size_t iii=0;iii<nshifts;++iii)
    {
        wsqPtrReal[iii]=wsq(iii).real();
        dzdtPtrReal[iii]=dzdt(iii).real();
        wsqPtrImag[iii]=wsq(iii).imag();
        dzdtPtrImag[iii]=dzdt(iii).imag();
    }
}