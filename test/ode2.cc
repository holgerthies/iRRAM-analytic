#include "iRRAM.h"
#include "IVPSOLVER.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;
using std::vector;

// F(t, y1, y2) = y2+1
REAL series1(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  if(v0 == 0 && v1 == 0 && v2 == 0) return 1;
  if(v0 == 0 && v1 == 0 && v2 == 1) return 1;
  return 0;
}

// F(t, y1, y2) = y2
REAL series1p(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  if(v0 == 0 && v1 == 0 && v2 == 1) return 1;
  return 0;
}
// F(t, y1, y2) = -y1
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
REAL FUN2p(const REAL& t){
  return cos(t);
}
void check(REAL x1, REAL y, int prec)
{
    iRRAM::cout << "result: " << endl;
    rwrite(x1, prec);
    iRRAM::cout << endl;
    iRRAM::cout << "should be " << endl;
    rwrite(y,prec);
    iRRAM::cout << endl;
    if(!bound(abs(x1-y), -prec)){
      iRRAM::cout << "ERRROR" << endl;
    } else {
      iRRAM::cout << "OK!" << endl;
    }
}

void solve_stepwise(const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f1, const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f2,const REAL& stepsize, const int steps, const int prec)
{
  auto f1c = f1;
  auto f2c = f2;
  auto r1=f1->to_analytic()->get_r();
  auto r2=f2->to_analytic()->get_r();
  REAL x=0, Y1=0, Y2=0;
  for(int i=1; i<=steps; i++){
    x += stepsize;
    auto F = ivp_solve({f1c,f2c});
    std::cout << "iteration "<<i<<endl;
    iRRAM::cout << "solving ode..." << endl;
    iRRAM::cout << "computing F1("<<x<<")..."<<endl;
    Y1 += F[0]->evaluate(stepsize);
    iRRAM::cout << "computing F2("<<x<<")..."<<endl;
    Y2 += F[1]->evaluate(stepsize);
    REAL new_r1 = r1-maximum(abs(x), maximum(abs(Y1), abs(Y2)));
    REAL new_r2 = r2-maximum(abs(x), maximum(abs(Y1), abs(Y2)));
    f1c = continuation(f1, 200000, new_r1,x,Y1,Y2);
    f2c = continuation(f2, 200000, new_r2,x,Y1,Y2);
    check(Y1, FUN1(x), prec);
    check(Y2, FUN2(x), prec);
    
  }
}

void compute(){
	
  auto f1 = make_analytic<REAL,REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1),4000000, 4000000);
  auto f1p = make_analytic<REAL,REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1p),4000000, 4000000);
  auto f1c = continuation(f1p, 200000, 200000,REAL(0),REAL(0),REAL(1));
  auto f2 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series2),4000000, 4000000);
  auto f2c = continuation(f2, 200000, 200000,REAL(0),REAL(0),REAL(1));
  int l,prec;
  //iRRAM::cin >>l>> prec;
  prec=50;
  l=100;
  // f continuation prec
  REAL x1= REAL(1)/REAL(16*l);
  REAL x2= REAL(1)/REAL(8*l);
  REAL x3= REAL(1)/REAL(10*l);

  iRRAM::cout << "computing f1("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  REAL y = f1c->evaluate(x1,x2,x3);
  REAL sol=fun1(x1,x2,x3);
  check(y, sol, prec);

  iRRAM::cout << "computing f2("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  y = f2c->evaluate(x1,x2,x3);
  sol=fun2(x1,x2,x3);
  check(y, sol, prec);
  
  solve_stepwise(f1,f2,0.1,100,prec);
  
  
  
}
