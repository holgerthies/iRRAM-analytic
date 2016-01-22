#include "iRRAM.h"
#include "ANALYTIC.h"
#include "ADDITION.h"
#include "coefficient_computation.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;

REAL sinseries1d(unsigned long n){
  if(n % 2 == 0) return 0;
  if((n-1) % 4 == 0) return inv_factorial(n);
  return -inv_factorial(n);
}
// sin(x)+1
REAL sinseriesmod1d(unsigned long n){
  if(n==0) return 1;
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

// series for sin(x1*x2*x3)
REAL sinseries2d(unsigned long n, unsigned long m){
  if(n != m ) return 0;
  if(n%2 == 0)
    return 0;
  else {
    if (0 == (n-1)%4)
      return inv_factorial(n);
    else
      return -inv_factorial(n);
  }
}
// series for sin(x1*x2*x3)+1
REAL sinseriesmod(unsigned long n, unsigned long m, unsigned long q){
  if(n+m+q==0) return 1;
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
  //auto sin_pwr = make_series<1,REAL>(std::function<REAL(const REAL&)>([] (const REAL& x) {return sin(x);}), 2,2);
	
  // auto sin_pwr3 = make_series<2,REAL>(std::function<REAL(const REAL&, const REAL&)>([] (const REAL& x, const REAL& y) {return sin(x*y);}), 2,2);
  //cout << get_coeff(*sin_pwr3, 1,1) << std::endl;
  auto f = make_analytic<REAL,REAL,REAL,REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(sinseries), 2, 2);
  auto g = make_analytic<REAL,REAL>(std::function<REAL(unsigned long)>(sinseries1d), 2,2);
  // ANALYTIC<1,REAL> h( 2,2, sin_pwr);
  // ANALYTIC<2,REAL> h3( 2,2,sin_pwr3);
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

  REAL y = g->evaluate(x1);
  REAL sol=sin(x1);
  iRRAM::cout << "checking sin(x)" << endl;
  checkResult(y, sol, prec);
  
  // y = h(x1);
  // sol=sin(x1);
  // iRRAM::cout << "checking sin(x) from function" << endl;
  // checkResult(y, sol, prec);

  y = f->evaluate(x1,x2,x3);
  sol=sin(x1*x2*x3);
  iRRAM::cout << "checking sin(x1*x2*x3)" << endl;
  checkResult(y, sol, prec);

  // // y = h3(x1,x2);
  // // sol=sin(x1*x2);
  // // iRRAM::cout << "checking sin(x1*x2) from function" << endl;
  // // chec

  // auto comp = compose(g,g);
  //  iRRAM::cout << "checking composition (1d)" << endl;
  //  y = comp(x1);
  //  sol=sin(sin(x1));
  //  checkResult(y, sol, prec);

  // auto comp2 = compose(g,f);
  // iRRAM::cout << "checking composition (3d)" << endl;
  // y = comp2(x1,x2,x3);
  // sol=sin(sin(x1*x2*x3));
  // checkResult(y, sol, prec);

  auto sum = f+f;
  y = sum->evaluate(x1,x2,x3);
  iRRAM::cout << "checking addition" << endl;
  sol=2*sin(x1*x2*x3);
  checkResult(y, sol, prec);

  auto sum2 = ( f+f )->to_analytic();
  y = sum2->evaluate(x1,x2,x3);
  iRRAM::cout << "checking addition by powerseries" << endl;
  sol=2*sin(x1*x2*x3);
  checkResult(y, sol, prec);
   
  //  auto prod = f*sum;
  //  y = prod(x1,x2,x3);
  //  iRRAM::cout << "checking multiplication" << endl;
  //  sol=(2*sin(x1*x2*x3)+2)*sin(x1*x2*x3);
  //  checkResult(y, sol, prec);
  //  ANALYTIC<1,REAL> gmod(std::function<REAL(unsigned long)>(sinseriesmod1d), 3,2);
  //  ANALYTIC<3,REAL> fmod(std::function<REAL(unsigned long, unsigned long, unsigned long)>(sinseriesmod), 3,2);


  // iRRAM::cout << "checking division (1d)" << endl;
  // // sin(x)/(sin(x)+1)
  // auto inv_test = g/gmod;   
  // y = inv_test(x1);
  // sol=sin(x1)/(1+sin(x1));
  // checkResult(y, sol, prec);


  // iRRAM::cout << "checking derivative (1d)" << endl;
  // // sin(x)
  // auto d1 = derive(g,0, 2);
  // y = d1(x1);
  // sol=-sin(x1);
  // checkResult(y, sol, prec);

  // iRRAM::cout << "checking derivative (3d)" << endl;
  // auto d2 = derive(derive(f,1, 1), 2, 2);
  // y = d2(x1, x2,x3);
  // sol=-cos(x1*x2*x3)*x1*x1*x1*x2*x2*x3-sin(x1*x2*x3)*x1*x1*x2-x1*x1*x2*sin(x1*x2*x3);
  // checkResult(y, sol, prec);


  // // check predefined functions
  // auto mon1=AnalyticFunction::monomial<3,2,REAL>(5); // x2^5
  // auto mon2=AnalyticFunction::monomial<3,1,REAL>(3); // x1^3
  // auto mon3=AnalyticFunction::monomial<3,3,REAL>(2); // x3^2
  // auto mon4=AnalyticFunction::monomial<3,1,REAL>(1); // x1

  // iRRAM::cout << "checking monomial (3d)" << endl;
  // y = mon1(x1, x2,x3);
  // sol=power(x2,5);
  // checkResult(y, sol, prec);
  // // x1^3*x3^2 - x2^5/(x1+1)
  // auto composed = mon2*mon3-mon1/(ANALYTIC<3,REAL>(1)+mon4);
  // y = composed(x1, x2,x3);
  // sol=power(x1,3)*power(x3,2)-power(x2,5)/(1+x1);
  // checkResult(y, sol, prec);

  // iRRAM::cout << "checking division (3d)" << endl;

  // // sin(x*y*z)/(sin(x*y*z)+1)
  // auto inv_test2 = f/fmod;   
  // y = inv_test2(x1,x2,x3);
  // sol=sin(x1*x2*x3)/(1+sin(x1*x2*x3));
  // checkResult(y, sol, prec);


}
