#include "iRRAM.h"
#include "IVPSOLVER.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;
using std::vector;

// f(x,y,z) = z
REAL series1(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  if(v0 == 0 && v1 == 0 && v2 == 0) return 1;
  if(v0 == 0 && v1 == 0 && v2 == 1) return 1;
  return 0;
}

// f(x,y,z) = y
REAL series2(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  if(v0 == 0 && v1 == 1 && v2 == 0) return -1;
  return 0;
}

REAL fun1(const REAL& t, const REAL& y1, const REAL& y2){
  return y2+1;
}

REAL fun2(const REAL& t, const REAL& y1, const REAL& y2){
  return -y1;
}

REAL FUN1(const REAL& t){
  return sin(t);
}

REAL FUN2(const REAL& t){
  return cos(t)-1;
}
void compute(){
	
  auto f1 = make_analytic<REAL,REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1),1, 2);
  auto f2 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series2),1, 2);
  auto F = ivp_solve({f1,f2});
  int l,prec;
  //iRRAM::cin >>l>> prec;
  prec=150;
  l=10;
  // f continuation prec
  REAL x1= REAL(1)/REAL(16*l);
  REAL x2= REAL(1)/REAL(8*l);
  REAL x3= REAL(1)/REAL(10*l);
  iRRAM::cout << "computing f1("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  REAL y = f1->evaluate(x1,x2,x3);
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  REAL sol=fun1(x1,x2,x3);
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
  iRRAM::cout << "computing f2("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  y = f2->evaluate(x1,x2,x3);
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  sol=fun2(x1,x2,x3);
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
  iRRAM::cout << "solving ode..." << endl;
  iRRAM::cout << "computing F1("<<x1<<")..."<<endl;
  REAL Y = F[0]->evaluate(x1);
  iRRAM::cout << "result: " << endl;
  rwrite(Y, prec);
  iRRAM::cout << endl;
  REAL SOL=FUN1(x1);
  iRRAM::cout << "should be " << endl;
  rwrite(SOL,prec);
  iRRAM::cout << endl;
  if(!bound(abs(SOL-Y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
  iRRAM::cout << "computing F2("<<x1<<")..."<<endl;
  Y = F[1]->evaluate(x1);
  iRRAM::cout << "result: " << endl;
  rwrite(Y, prec);
  iRRAM::cout << endl;
  SOL=FUN2(x1);
  iRRAM::cout << "should be " << endl;
  rwrite(SOL,prec);
  iRRAM::cout << endl;
  if(!bound(abs(SOL-Y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
}
