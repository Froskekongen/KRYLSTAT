#define EIGEN_YES_I_KNOW_SPARSE_MODULE_IS_NOT_STABLE_YET

#include <iostream>
#include "eigen_det_app.h"
#include "../../eigen_pow.h"

using namespace ColPack;
int main()
{
  
  int_type sizeAmat=20;
  int_type nshifts=12;
  // int_type Nsamples=10;
  
  std::cout << "We have come to line 15 \n\n";
  
  SparseMatrix<double,RowMajor,int> AB(sizeAmat,sizeAmat);
  AB.reserve(3*sizeAmat);
  for (unsigned int iii=0;iii<sizeAmat;++iii){AB.insert(iii,iii)=2.1;}
  for (unsigned int iii=0;iii<sizeAmat-1;++iii){AB.insert(iii+1,iii)=-1.;}
  for (unsigned int iii=0;iii<sizeAmat-1;++iii){AB.insert(iii,iii+1)=-1.;}
  AB.finalize();
  SparseMatrix<double,RowMajor,int> A(AB);
  SparseTriangularView< SparseMatrix<double,RowMajor,int>, Upper > CCC(A);
  
  std::cout << CCC;
  
  // std::cout << A;
  
  // MappedSparseMatr
  // Matrix<double,Dynamic,Dynamic> samples(sizeAmat,Nsamples);
  
  //SparseMatrix<double,RowMajor,unsigned int> BOLLE(A*A); //This is a problem line!
  std::cout << "We have come to line 27 \n\n";
  
  
  
  
  
  double minEig=0.01,maxEig=5.;
  double cg_tol=1e-4;
  int_type cg_maxiter=2000;
  double logdet=0;
  int_type kdist=5;
  double logdet2;
  
   Matrix<double , Dynamic,1> BB(static_cast<std::size_t>(1));
//   
//   unsigned int *BOB1=A._innerIndexPtr();
//   unsigned int *BOB2=A._outerIndexPtr();
  
  //SpStructCont<unsigned int> SpS(BOB2,BOB1,sizeAmat);
  
   logDetApp(A, minEig, maxEig, nshifts, cg_tol,cg_maxiter, kdist, logdet,BB,3);
   logDetMargVar(A, cg_tol, cg_maxiter, kdist, logdet2, BB,3);
  
  std::cout << "Logdet2: " << logdet2 << "\n\n";
  std::cout << "Logdet: " << logdet << "\n\n";
  
  // MatrixXd MM(A);
  // std::cout << MM;
  
  
  
//   double *intConst = new double;
//   std::complex<double> *wsq = new std::complex<double>[nshifts];
//   std::complex<double> *dzdt = new std::complex<double>[nshifts];
//   
//   logIntPoints(nshifts, minEig, maxEig, intConst, wsq, dzdt);
//   
//   for (int iii=0;iii<nshifts;++iii)
//   {
//     std::cout << "wsq(" << iii<< ")=" << wsq[iii] << "        dzdt(" << iii << ")=" << dzdt[iii] << ")       intConst=" << *intConst << "\n\n";
//   }

//   {
//     int_type powwer=3;
//     SparseMatrix<double> B=A;
//     //B.reserve(6*sizeAmat);
//     eigen_power_sparse(A,powwer,B);
//     
//     int *col_indices=B._innerIndexPtr();
//     int *row_offsets=B._outerIndexPtr();
//   
//     SpStructCont<int> SPS(row_offsets,col_indices,sizeAmat);
//   
//     SPS.printSPS();
//     GraphColoringInterface * ColOr = new GraphColoringInterface(SRC_MEM_ADOLC, SPS.pptr, SPS.num_rows);
//     ColOr->Coloring("NATURAL", "DISTANCE_TWO");
//     ColOr->PrintVertexColoringMetrics();
//     delete ColOr;
//   }
  
  VectorXd z=VectorXd::Ones(sizeAmat);
  DiagonalWrapper<const VectorXd> M=z.asDiagonal();
  
  VectorXd margVars=VectorXd::Zero(sizeAmat);
  VectorXi probes=VectorXi::Zero(sizeAmat);
  
  computeProbingVectors(A, probes, 1, 10,3);
  
  
  computeMargVars(A, M, probes, margVars, cg_tol, cg_maxiter);
  
  std::cout << "\n\n" << margVars << "\n\n";
  
  
  
  
  

  return 1;
}