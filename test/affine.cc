#include "iRRAM.h"
#include "AAREAL.h"
#include "TaylorModel.h"
#include <iostream>
using namespace std;
using namespace iRRAM;
REAL to_real(const AAREAL& x)
{
  return x.to_real();
}
REAL to_real(const REAL& x)
{
  return x;
}
// typedef REAL RTYPE;
typedef AAREAL RTYPE;
void compute(){
  int method;
  static int it=0;
  it++;
  std::cout << "Iter "<<it<<std::endl;
  
  bool compare=false;
  iRRAM::cin >> method;
  const int steps = 150000;
  const int prec=1000;
  REAL a=3.6, x0="0.4";
  if(method == 0){
    AAREAL aa(a);
    AAREAL x(x0);
    REAL xr(x0);
    AAREAL one(1);
    
    for(int i=0; i<steps; i++){
      x = a*x*(one-x);
      x.clean();
      if(compare)
        xr = a*xr*(1-xr);
      if(compare || i == steps-1)
      {
        iRRAM::cout << "iter " << (i+1) << ": ";
      REAL r=to_real(x);
      // x.print_error();
      
      iRRAM::cout << std::endl;
      rwrite(r, prec);
       iRRAM::cout << std::endl;
      if(compare)
      {
        rwrite(xr, prec);
      iRRAM::cout << std::endl;
      }
      sizetype error;
      r.geterror(error);
      std::cout << "SUM AA: " << error.mantissa << ", " << error.exponent << std::endl;
      if(compare)
      {
        xr.geterror(error);
       std::cout << "SUM R: " << error.mantissa << ", " << error.exponent << std::endl;
      }
      }
      
    }
  }
  else if(method == 1){
    TM aa(a);
    TM x(x0);
    REAL xr(x0);
    TM one(1);
    
    for(int i=0; i<steps; i++){
      x = aa*x*(one-x);
      if(compare)
        xr = a*xr*(1-xr);
      if(i % 2 == 0) x.round();
      
      if(compare || i == steps-1)
      {
        iRRAM::cout << "iter " << (i+1) << ": ";
      REAL r=REAL(x);
      // x.print_error();
      
      iRRAM::cout << std::endl;
      rwrite(r, prec);
       iRRAM::cout << std::endl;
      if(compare)
      {
        rwrite(xr, prec);
      iRRAM::cout << std::endl;
      }
      sizetype error;
      r.geterror(error);
      std::cout << "SUM AA: " << error.mantissa << ", " << error.exponent << std::endl;
      if(compare)
      {
        xr.geterror(error);
       std::cout << "SUM R: " << error.mantissa << ", " << error.exponent << std::endl;
      }
      }
      
    }
}
  else
  {
    REAL aa(a);
    REAL x(x0);
    REAL one(1);
    for(int i=0; i<steps; i++){
      x = aa*x*(one-x);
    }
    iRRAM::cout << "iter " << steps << ": ";
    iRRAM::cout << std::endl;
    rwrite(x, prec);
    iRRAM::cout << std::endl;
    sizetype error;
    x.geterror(error);
    std::cout << "SUM: " << error.mantissa << ", " << error.exponent << std::endl;
  }
  
  
}

