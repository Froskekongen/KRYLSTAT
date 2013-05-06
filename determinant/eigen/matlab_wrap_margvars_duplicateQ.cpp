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
#include "eigen_det_app.h"
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

using namespace std;

typedef MappedSparseMatrix<double,RowMajor,size_t> MappedSparseMatrixType;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
    {
        mexErrMsgTxt("Need 4 input arguments: matrix (upper triangular), probingvector, cgTol, cgMaxIter");
    }
    if (nlhs>1)
    {
        mexErrMsgTxt("Maximum of 2 output arguments");
    }
    
    int i=0;
    MappedSparseMatrixType B(mxGetM(prhs[i]), mxGetM(prhs[i]), mxGetNzmax(prhs[i]), mxGetJc(prhs[i]), mxGetIr(prhs[i]), mxGetPr(prhs[i]));
    SparseSelfAdjointView< MappedSparseMatrixType , Lower> A(B); //Lower corresponds to triu (upper) in matlab
    
    
    int *probePointer;
    probePointer = (int*)mxGetData(prhs[1]);
    Map<VectorXi> eigenProbe(probePointer,A.rows());
    
    
    double cgTol=mxGetScalar(prhs[2]);
    size_t cgMaxIter(static_cast<int>(mxGetScalar(prhs[3])));
    
    plhs[0]=mxCreateDoubleMatrix(A.cols(),1,mxREAL);
    double *dataPointer;
    dataPointer=mxGetPr(plhs[0]);
    Map<VectorXd> probeAinv(dataPointer,A.cols());
    
    computeMargVarsDuplicateQ(A,eigenProbe,probeAinv,cgTol,cgMaxIter);
    
    
    
    
}