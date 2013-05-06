#include <stdio.h>
#include <iostream>
#include "sncndn.h"

int main()
{
  
  
  double sn,cn,dn;
  double u(0.2);
  double m(0.99);
  
  recursiveSNCNDN(u,m,sn,cn,dn);
  
  double K,Kp;
  ellKKP(u,K,Kp);
  
  /*std::cout << "Real parts: sn=" << sn.real << ", cn=" << cn.real << ", dn=" << dn.real << "\n";
  std::cout << "Imaginary parts: sn=" << sn.imag << ", cn=" << cn.imag << ", dn=" << dn.imag << "\n";*/
  
  std::cout << "sn=" << sn << ", cn=" << cn << ", dn=" << dn << "\n";
  std::cout << "K=" << K << ", Kp=" << Kp << "\n";
  
  std::cout << mp_real::_eps << "\n \n";
  
  
  const int N=7;
  double *intConst=new double;
  cpxDbl *wsq=new cpxDbl[N];
  cpxDbl *dzdt=new cpxDbl[N];
  
  logIntPoints(N,0.0001,20.001,intConst,wsq,dzdt);
  
  std::cout << "IntConst=" << *intConst << "\n";
  for (int iii=0;iii<N;++iii)
  {
    std::cout << "wsq(" << iii << ")=" << wsq[iii] << ", dzdt(" << iii <<")=" << dzdt[iii] << "\n";
  }
  
  double *wsqSq=new double[N];
  double *dzdtSq=new double[N];
  sqrtIntPoints(N,0.0001,20.001,intConst,wsqSq,dzdtSq);
  
  
  std::cout << "IntConst=" << *intConst << "\n";
  for (int iii=0;iii<N;++iii)
  {
    std::cout << "wsq(" << iii << ")=" << wsqSq[iii] << ", dzdt(" << iii <<")=" << dzdtSq[iii] << "\n";
  }
  
  
  
  
  delete intConst;
  delete[] wsq;
  delete[] dzdt;
  delete[] wsqSq;
  delete[] dzdtSq;
  
  
  return 1;
}