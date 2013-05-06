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
#include "eigen_det_app.h"
#include <math.h>
#include "mex.h"
#include "matrix.h"
#include <stdlib.h>

using namespace std;

typedef MappedSparseMatrix<double,RowMajor,size_t> MappedSparseMatrixType;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 8)
    {
        mexErrMsgTxt("Need 8 input arguments: matrix (upper triangular), probingvector, wsq, dzdt, intconst, cgTol, cgMaxIter, verbosity");
    }
    if (nlhs>2)
    {
        mexErrMsgTxt("Maximum of 2 output arguments");
    }
    
     double ldet=0;
     
     int i=0;
    MappedSparseMatrixType B(mxGetM(prhs[i]), mxGetM(prhs[i]), mxGetNzmax(prhs[i]), mxGetJc(prhs[i]), mxGetIr(prhs[i]), mxGetPr(prhs[i]));
    //SparseMatrix<double,RowMajor,std::size_t> C(B);
    SparseSelfAdjointView< MappedSparseMatrixType , Lower> A(B); //Lower corresponds to triu (upper) in matlab
  
    
    int *probePointer;
    probePointer = (int*)mxGetData(prhs[1]);
    Map<VectorXi> eigenProbe(probePointer,A.rows());
    
    double *wsqPtrReal=mxGetPr(prhs[2]), *wsqPtrImag=mxGetPi(prhs[2]);
    double *dzdtPtrReal=mxGetPr(prhs[3]), *dzdtPtrImag=mxGetPi(prhs[3]);
    double intConst=mxGetScalar(prhs[4]); 
    double cgTol=mxGetScalar(prhs[5]);
    size_t cgMaxIter(static_cast<int>(mxGetScalar(prhs[6])));
    int verbosity(mxGetScalar(prhs[7]));
     
    int nshifts=mxGetM(prhs[2]);
    
    VectorXcd wsqEig(nshifts),dzdtEig(nshifts);
    wsqEig.real() = Map<VectorXd>(wsqPtrReal,nshifts);
    wsqEig.imag() = Map<VectorXd>(wsqPtrImag,nshifts);
    dzdtEig.real() = Map<VectorXd>(dzdtPtrReal,nshifts);
    dzdtEig.imag() = Map<VectorXd>(dzdtPtrImag,nshifts);
    
    VectorXd mv(1);
    
    // std::cout << wsqEig << "\n"<< dzdtEig << "\n" << intConst << "\n" << eigenProbe << "\n" << cgTol << "\n" << cgMaxIter << "\n\n";
    
    logDetBase(A,eigenProbe,wsqEig,dzdtEig,intConst,cgTol,cgMaxIter,ldet,mv,verbosity);   
    
    
   plhs[0]=mxCreateDoubleScalar(ldet);
}