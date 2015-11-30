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


REAL sinseries1d(unsigned long n){
  if(n % 2 == 0) return 0;
  if((n-1) % 4 == 0) return inv_factorial(n);
  return -inv_factorial(n);
}

// series for sin(x1*x2*x3)
REAL sinseries(unsigned long n, unsigned long m, unsigned long q){
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

void checkResult(const REAL& y,const REAL& sol, const int prec ){
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
}

void compute(){
	
  ANALYTIC<3,REAL> f(std::function<REAL(unsigned long, unsigned long, unsigned long)>(sinseries), 2, 2);
  ANALYTIC<1,REAL> g(std::function<REAL(unsigned long)>(sinseries1d), 2,2);
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

  REAL y = g(x1);
  REAL sol=sin(x1);
  iRRAM::cout << "checking sin(x)" << endl;
  checkResult(y, sol, prec);

  y = f(x1,x2,x3);
  sol=sin(x1*x2*x3);
  iRRAM::cout << "checking sin(x1*x2*x3)" << endl;
  checkResult(y, sol, prec);

  auto sum = f+f+ANALYTIC<3,REAL>(2);
  y = sum(x1,x2,x3);
  iRRAM::cout << "checking addition" << endl;
  sol=2*sin(x1*x2*x3)+2;
  checkResult(y, sol, prec);

  auto prod = f*sum;
  y = prod(x1,x2,x3);
  iRRAM::cout << "checking multiplication" << endl;
  sol=(2*sin(x1*x2*x3)+2)*sin(x1*x2*x3);
  checkResult(y, sol, prec);


  auto comp = compose(g,g);
  iRRAM::cout << "checking composition (1d)" << endl;
  y = comp(x1);
  sol=sin(sin(x1));
  checkResult(y, sol, prec);

  auto comp2 = compose(g,f);
  iRRAM::cout << "checking composition (3d)" << endl;
  y = comp2(x1,x2,x3);
  sol=sin(sin(x1*x2*x3));
  checkResult(y, sol, prec);
}
