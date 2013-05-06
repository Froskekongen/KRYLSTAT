#include <stdio.h>
#include <iostream>
#include <complex>
#include "jacobiEllipticDouble.h"


int main()
{
  std::complex<double> sn,cn,dn, u(0.2),m(0.99);
  double dummy;
  
  jacobiEllipticDouble::recursiveSNCNDN(u,m,sn,cn,dn,dummy);
  
  double K,Kp;
  jacobiEllipticDouble::ellKKP(u.real(),K,Kp);
  
  /*std::cout << "Real parts: sn=" << sn.real() << ", cn=" << cn.real() << ", dn=" << dn.real << "\n";
  std::cout << "Imaginary parts: sn=" << sn.imag() << ", cn=" << cn.imag() << ", dn=" << dn.imag() << "\n";*/
  
  std::cout << "sn=" << sn << ", cn=" << cn << ", dn=" << dn << "\n";
  std::cout << "K=" << K << ", Kp=" << Kp << "\n";
  
  // td::cout << mp_real::_eps << "\n \n";
  
  
  const int N=7;
  double intConst(0);
  VectorXcd wsq=VectorXcd::Zero(N);
  VectorXcd dzdt=VectorXcd::Zero(N);
  
  jacobiEllipticDouble::logIntPoints(N,0.0001,20.001,intConst,wsq,dzdt,3);
  
//   std::cout << "IntConst=" << intConst << "\n";
//   for (int iii=0;iii<N;++iii)
//   {
//     std::cout <<"wsq(" << iii << ")=" << wsq(iii) << ", dzdt(" << iii <<")=" << dzdt(iii) << "\n";
//   }
  
  double intConstSq(0);
  VectorXd wsqSq=VectorXd::Zero(N);
  VectorXd dzdtSq=VectorXd::Zero(N);
  jacobiEllipticDouble::sqrtIntPoints(N,0.0001,20.001,intConstSq,wsqSq,dzdtSq,3);
  
  return 1;
}