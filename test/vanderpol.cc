#include "iRRAM.h"
#include "ANALYTIC.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;
using std::vector;
REAL get_mu()
{
  return 10;
  
}
REAL series1(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  if(v0 == 0 && v1 == 0 && v2 == 0) return 0.1;
  if(v0 == 0 && v1 == 0 && v2 == 1) return 1;
  return 0;
}

REAL series2(const unsigned long v0, const unsigned long v1, const unsigned long v2){
  const REAL mu=get_mu();
  if(v0 == 0 && v1 == 0 && v2 == 0) return 0.1*mu;
  if(v0 == 0 && v1 == 0 && v2 == 1) return mu;
  if(v0 == 0 && v1 == 1 && v2 == 0) return -1;
  if(v0 == 0 && v1 == 2 && v2 == 0) return -0.1*mu;
  if(v0 == 0 && v1 == 2 && v2 == 1) return -mu;
  return 0;
}

REAL fun1(const REAL& t, const REAL& y1, const REAL& y2){
  return y2+REAL(0.1);
}

REAL fun2(const REAL& t, const REAL& y1, const REAL& y2){
  const REAL mu=get_mu();
  return mu*(1-y1*y1)*(y2+REAL(0.1))-y1;
}

void solve_stepwise(const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f1, const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f2,const REAL& stepsize, const int steps, const int prec)
{
  auto f1c = f1;
  auto f2c = f2;
  auto r1=f1->to_analytic()->get_r();
  auto r2=f2->to_analytic()->get_r();
  auto M1 = f1->to_analytic()->get_M();
  auto M2 = f2->to_analytic()->get_M();
  REAL x=0, Y1=0, Y2=0.1;
  for(int i=1; i<=steps; i++){
    x += stepsize;
    auto F = ivp_solve({f1c,f2c});
    Y1 += F[0]->evaluate(stepsize);
    Y2 += F[1]->evaluate(stepsize);
    REAL new_r1 = r1-maximum(abs(x), maximum(abs(Y1), abs(Y2)));
    REAL new_r2 = r2-maximum(abs(x), maximum(abs(Y1), abs(Y2)));
    f1c = continuation(f1, M1, new_r1,x,Y1,Y2);
    f2c = continuation(f2, M2, new_r2,x,Y1,Y2);
    iRRAM::cout << x << " " << Y1 << " " << Y2 << std::endl;
    
  }
}

void compute(){
	
  REAL r = 1;
  REAL mu = get_mu();
  REAL M1 = r+REAL(0.1);
  REAL M2= mu*(abs(r*r-1))*(r+REAL(0.1))+r;
  auto f1 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1),M1, r);
  auto f2 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series2),M2, r);
  int prec=2;
  solve_stepwise(f1, f2, 0.1*r/M2,100, prec);
}
