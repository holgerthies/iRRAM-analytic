#include "iRRAM.h"
#include "ANALYTIC.h"
using namespace iRRAM;
using std::endl;
REAL inv_factorial(const int n){
  if ((n!=0)&&(n*log(n)-n > 2*-ACTUAL_STACK.actual_prec)){
    REAL return_value(0);
    sizetype error;
    sizetype_set(error,1,ACTUAL_STACK.actual_prec);
    return_value.seterror(error);
    return return_value;
  }
  if (n==0)
    return REAL(1);
  REAL inv_fact=inv_factorial(n-1)/REAL(n);
  return inv_fact;
}


// sin(2*pi*x)
REAL sinseries(const std::vector<unsigned long>& v){
  int n=v[0];
  int m=v[1];
  int q=v[2];
  if(n != m || n != q) return 0;
  if(n%2 == 0)
    return 0;
  else {
    if (0 == (n-1)%4)
      return inv_factorial(n);
    else
      return -inv_factorial(n);
  }
}

void compute(){
	
  ANALYTIC<3,REAL> f(std::shared_ptr<std::function<REAL(const std::vector<unsigned long>&)>>(new std::function<REAL(const std::vector<unsigned long>&)>(sinseries)), 2, 2);

  //auto g = AnalyticFunction::projection<3,1,REAL>();
  //g = power(g,5);
  int l,prec;
  iRRAM::cin >>l>> prec;
  // rwrite(g({0.2,0.3, 0.5}), prec);
  iRRAM::cout << endl;
  // f continuation prec
  REAL x1= REAL(1)/REAL(4*l);
  REAL x2= REAL(1)/REAL(2*l);
  REAL x3= REAL(1)/REAL(8*l);
  REAL y = f({x1,x2,x3});
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  REAL sol=sin(x1*x2*x3);
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
  auto sum = f+f;
  y = sum({x1,x2,x3});
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  sol=2*sin(x1*x2*x3);
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
  auto prod = f*sum;
  y = prod({x1,x2,x3});
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  sol=2*sin(x1*x2*x3)*sin(x1*x2*x3);
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
}
