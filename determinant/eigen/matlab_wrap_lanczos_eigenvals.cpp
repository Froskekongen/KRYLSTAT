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

#include <Eigen/Eigen>
#include <iostream>
#include "eigen_det_app.h"
#include "../../lanczos/eigen/eigen_lanczos.h"
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
        mexErrMsgTxt("Need 5 input arguments: matrix (upper triangular), cgTol, cgMaxIter, verbosity");
    }
    if (nlhs>2)
    {
        mexErrMsgTxt("Maximum of 2 output arguments");
    }
    
    int i=0;
    MappedSparseMatrixType B(mxGetM(prhs[i]), mxGetM(prhs[i]), mxGetNzmax(prhs[i]), mxGetJc(prhs[i]), mxGetIr(prhs[i]), mxGetPr(prhs[i]));
    //SparseMatrix<double,RowMajor,std::size_t> C(B);
    SparseSelfAdjointView< MappedSparseMatrixType , Lower> A(B); //Lower corresponds to triu (upper) in matlab
    
    double cgTol=mxGetScalar(prhs[1]);
    size_t cgMaxIter(static_cast<int>(mxGetScalar(prhs[2])));
    int verbosity(mxGetScalar(prhs[3]));
    
    double minEig(0), maxEig(0);
    
    
    VectorXd b=VectorXd::Random(A.cols());
    VectorXd x=VectorXd::Zero(A.cols());
    default_monitor<double> mon(b,cgMaxIter,cgTol);
    
    // verbosity currently not implemented
    eigen_lanc_extremal_eig_mult(A,b,x,minEig,maxEig,mon);
    
    plhs[0]=mxCreateDoubleScalar(minEig);
    plhs[1]=mxCreateDoubleScalar(maxEig);
        
}
    