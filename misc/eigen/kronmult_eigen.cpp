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

#include <iostream>
#include "kronMatVec.h"
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

using namespace Eigen;
using namespace std;

typedef MappedSparseMatrix<double,ColMajor,size_t> MappedSparseMatrixType;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 2)
    {
        mexErrMsgTxt("Need 2 input arguments: cell array of matrices, inputVector");
    }
    if (nlhs>1)
    {
        mexErrMsgTxt("Maximum of 1 output arguments");
    }
    int i=0;
    int nbs = mxGetNumberOfElements (prhs[i]);
    std::vector<MappedSparseMatrixType, aligned_allocator<MappedSparseMatrixType> > Qs;
    for(int j = 0; j < nbs; j ++) 
    {
        mxArray *mwP = mxGetCell (prhs[i], j);
        Qs.push_back( MappedSparseMatrixType(  mxGetN(mwP), mxGetM(mwP), mxGetNzmax(mwP), mxGetJc(mwP), mxGetIr(mwP), mxGetPr(mwP) )  );
    }
    
    double *inpVec = mxGetPr(prhs[1]);
    int length=(int)mxGetM(prhs[1]);
    Map<VectorXd> inpEig(inpVec,length);
    
    // std::cout << inpEig[0] <<std::endl;
    
    plhs[0]=mxCreateDoubleMatrix(length,1,mxREAL);
    double *outVec = mxGetPr(plhs[0]);
    Map<VectorXd> outEig(outVec,length);
    
    
    krylstat_misc::kronMatVec(Qs,inpEig,outEig);
    
    
    
}