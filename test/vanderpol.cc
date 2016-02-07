#include "iRRAM.h"
#include "ANALYTIC.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;
using std::vector;
REAL series1(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  if(v0 == 0 && v1 == 0 && v2 == 1) return 1;
  return 0;
}

REAL series2(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  const REAL mu=0.1;
  if(v0 == 0 && v1 == 0 && v2 == 0) return -2;
  if(v0 == 0 && v1 == 0 && v2 == 1) return -3*mu;
  if(v0 == 0 && v1 == 1 && v2 == 0) return -1;
  if(v0 == 0 && v1 == 1 && v2 == 1) return -4*mu;
  if(v0 == 0 && v1 == 2 && v2 == 1) return -mu;
  return 0;
}

REAL fun1(const REAL& t, const REAL& y1, const REAL& y2){
  return y2;
}

REAL fun2(const REAL& t, const REAL& y1, const REAL& y2){
  const REAL mu=10;
  return mu*(1-(y1+2)*(y1+2))*y2-y1-2;
}

void compute(){
	
  auto f1 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1),1, 100);
  auto f2 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series2),1, 100);
  int prec=20;
  int l=100;
  // f continuation prec
  // REAL x1= REAL(1)/REAL(16*l);
  // REAL x2= REAL(1)/REAL(8*l);
  // REAL x3= REAL(1)/REAL(10*l);
  // iRRAM::cout << "computing f1("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  // REAL y = f1->evaluate(x1,x2,x3);
  // iRRAM::cout << "result: " << endl;
  // rwrite(y, prec);
  // iRRAM::cout << endl;
  // REAL sol=fun1(x1,x2,x3);
  // iRRAM::cout << "should be " << endl;
  // rwrite(sol,prec);
  // iRRAM::cout << endl;
  // if(!bound(abs(sol-y), -prec)){
  //   iRRAM::cout << "ERRROR" << endl;
  // } else {
  //   iRRAM::cout << "OK!" << endl;
  // }
  // iRRAM::cout << "computing f2("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  // y = f2->evaluate(x1,x2,x3);
  // iRRAM::cout << "result: " << endl;
  // rwrite(y, prec);
  // iRRAM::cout << endl;
  // sol=fun2(x1,x2,x3);
  // iRRAM::cout << "should be " << endl;
  // rwrite(sol,prec);
  // iRRAM::cout << endl;
  // if(!bound(abs(sol-y), -prec)){
  //   iRRAM::cout << "ERRROR" << endl;
  // } else {
  //   iRRAM::cout << "OK!" << endl;
  // }
  // iRRAM::cout << "solving ode..." << endl;
  auto F = ivp_solve({ f1,f2 });
  for(int i=0; i<=3*l;i++)
  {
    REAL x1=REAL(i)/REAL(l);
    rwrite(x1, prec);
    iRRAM::cout << " ";
    REAL Y = F[0]->evaluate(x1)+2;
    rwrite(Y, prec);
    iRRAM::cout << " ";
    Y = F[1]->evaluate(x1);
    rwrite(Y, prec);
    iRRAM::cout << endl;
  }
}
