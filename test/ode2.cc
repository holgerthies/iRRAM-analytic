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

// F(t, y1) = -y1
REAL series3(const unsigned long v0, const unsigned long v1){
  if(v0 == 0 && v1 == 1) return -1;
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
void check(REAL x1, REAL y10, int prec)
{
    iRRAM::cout << "result: " << endl;
    rwrite(x1, prec);
    iRRAM::cout << endl;
    iRRAM::cout << "should be " << endl;
    rwrite(y10,prec);
    iRRAM::cout << endl;
    if(!bound(abs(x1-y10), -prec)){
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

void euler_step(const std::shared_ptr<Node<REAL, REAL, REAL,REAL>>& f1,const std::shared_ptr<Node<REAL, REAL, REAL,REAL>>& f2, const REAL& t0, const REAL& y00,const REAL& y01, REAL& h, REAL& y10, REAL& y11)
{
  REAL M1=2*f1->to_analytic()->get_M()/f1->to_analytic()->get_r();
  REAL M2=2*f2->to_analytic()->get_M()/f2->to_analytic()->get_r();
  REAL M = maximum(M1, M2);
  h = 0.5;
  
  y10 = h*f1->evaluate(t0,y00, y01);
  y11 = h*f2->evaluate(t0,y00, y01);
  sizetype eval_error1, eval_error2, eval_error, trunc_error;
  y10.geterror(eval_error1);
  y11.geterror(eval_error2);
  sizetype_max(eval_error, eval_error1, eval_error2);
  REAL error = (h*h*M)/2;
  trunc_error = real_to_error(error);
  while (sizetype_less(eval_error, trunc_error) &&
             (trunc_error.exponent >= ACTUAL_STACK.actual_prec) ){
    h /= 2;
    y10 /= 2;
    y11 /= 2;
    y10.geterror(eval_error1);
    y11.geterror(eval_error2);
    sizetype_max(eval_error, eval_error1, eval_error2);
    error /= 4;
    trunc_error = real_to_error(error);
  }
  sizetype total_error1, total_error2;
  
  y10 += y00;
  y11 += y01;
  y10.geterror(eval_error1);
  y11.geterror(eval_error2);
  sizetype_add(total_error1, eval_error1, trunc_error);
  sizetype_add(total_error2, eval_error2, trunc_error);
  y10.seterror(total_error1);
  y11.seterror(total_error2);
}

void heun_step(const std::shared_ptr<Node<REAL, REAL, REAL,REAL>>& f1,const std::shared_ptr<Node<REAL, REAL, REAL,REAL>>& f2, const REAL& t0, const AAREAL& y00,const AAREAL& y01,const REAL& h, AAREAL& y10, AAREAL& y11)
{
  auto f1a = f1->to_analytic();
  auto f2a = f2->to_analytic();
  REAL M1=f1a->get_M(), r1=f1a->get_r(), M2=f2a->get_M(), r2 = f2a->get_r();
  REAL B1 = 2*power(M1,3)/(r1*r1);
  REAL B2 = 2*power(M2,3)/(r2*r2);
  REAL M = maximum(B1, B2);
  //h = 0.5;
  
  
  sizetype eval_error1, eval_error2, eval_error, trunc_error;
  REAL error;
  
  //do{
    //h /= 2;
    REAL ev1 = f1->evaluate(t0,y00.to_real(), y01.to_real());
    AAREAL y00p = y00+h*ev1;
    AAREAL y01p = y01+h*ev1;
    y10 = y00+(h/2)*(ev1+f1->evaluate(t0+h, y00p.to_real(), y01p.to_real()));
    REAL ev2 = f2->evaluate(t0,y00.to_real(), y01.to_real());
    AAREAL y10p = y00+h*ev2;
    AAREAL y11p = y01+h*ev2;
    y11 = y01+(h/2)*(ev2+f2->evaluate(t0+h, y10p.to_real(), y11p.to_real()));
    y10.to_real().geterror(eval_error1);
    y11.to_real().geterror(eval_error2);
    sizetype_max(eval_error, eval_error1, eval_error2);
    //error = power(h,3)*M;
    trunc_error = real_to_error(error);
  // } while (
  //          (trunc_error.exponent >= ACTUAL_STACK.actual_prec) );
  sizetype total_error1, total_error2;
  //y10.add_error(error);
  //y11.add_error(error);
  
  //y10.geterror(eval_error1);
  //y11.geterror(eval_error2);
  // sizetype_add(total_error1, eval_error1, trunc_error);
  // sizetype_add(total_error2, eval_error2, trunc_error);
  // y10.seterror(total_error1);
  // y11.seterror(total_error2);
}
void solve_stepwise_euler(const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f1, const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f2, const REAL& last, const int prec)
{
  REAL x=0, Y1=0, Y2=0;
  int i=0;
  while(x<last){
    REAL h, nY1, nY2;
    euler_step(f1,f2, x, Y1, Y2,h, nY1, nY2 );
    x += h;
    Y1 = nY1;
    Y2 = nY2;
    i++;
    if(i % 10000 == 0){
    std::cout << "iteration "<<i<<endl;
    iRRAM::cout << " Approx iterations needed " << (last/h) << endl;
    iRRAM::cout << "solving ode..." << endl;
    iRRAM::cout << "computing F1("<<x<<")..."<<endl;
    iRRAM::cout << "computing F2("<<x<<")..."<<endl;
    check(Y1, FUN1(x), prec);
    check(Y2, FUN2(x), prec);
    }
    
  }
}

REAL solve_stepwise_heun(const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f1, const std::shared_ptr<Node<REAL, REAL, REAL, REAL>>& f2, const REAL& last, const int prec)
{
  REAL x=0;
  AAREAL Y1, Y2;
  int i=0;
  std::cout << prec << std::endl;
  std::cout << ACTUAL_STACK.actual_prec << std::endl;
  
  while(x<last){
    REAL h=power(2, prec/3+1);
    AAREAL nY1, nY2;
    heun_step(f1,f2, x, Y1, Y2,h, nY1, nY2 );
    x += h;
    Y1 = nY1;
    Y2 = nY2;
    i++;
    if(i % 5 == 0){
      Y1.clean();
      Y2.clean();
    }
    // if(i % 10000 == 0 || true){
    // std::cout << "iteration "<<i<<endl;
    // iRRAM::cout << " Approx iterations needed " << (last/h) << endl;
    // iRRAM::cout << "solving ode..." << endl;
    // iRRAM::cout << "computing F1("<<x<<")..."<<endl;
    // iRRAM::cout << "computing F2("<<x<<")..."<<endl;
    // check(Y1.to_real(), FUN1(x), prec);
    // check(Y2.to_real(), FUN2(x), prec);
    // }
    
  }
  return Y1.to_real();
}

REAL limit_test(const int prec, const REAL& x)
{
  auto f1 = make_analytic<REAL,REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1),4, 4);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(f1)->add_algorithm(fun1);
  auto f2 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series2),4, 4);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(f2)->add_algorithm(fun2);
  return solve_stepwise_heun(f1, f2, x, prec);
}

void compute(){
	
  auto f1 = make_analytic<REAL,REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1),4, 4);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(f1)->add_algorithm(fun1);
  auto f1p = make_analytic<REAL,REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series1p),4, 4);
  auto f1c = continuation(f1p, 200000, 200000,REAL(0),REAL(0),REAL(1));
  auto f2 = make_analytic<REAL, REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(series2),4, 4);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(f2)->add_algorithm(fun2);
  auto f2c = continuation(f2, 200000, 200000,REAL(0),REAL(0),REAL(1));
  auto f3 = make_analytic<REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long)>(series3),2, 2);
  int l,prec;
  //iRRAM::cin >>l>> prec;
  prec=15;
  l=100;
  // f1 continuation prec
  REAL x1= REAL(1)/REAL(16*l);
  REAL x2= REAL(1)/REAL(8*l);
  REAL x3= REAL(1)/REAL(10*l);

  // iRRAM::cout << "computing f1("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  // REAL y10 = f1c->evaluate(x1,x2,x3);
  // REAL sol=fun1(x1,x2,x3);
  // check(y10, sol, prec);

  // iRRAM::cout << "computing f2("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
  // y10 = f2c->evaluate(x1,x2,x3);
  // sol=fun2(x1,x2,x3);
  // check(y10, sol, prec);

  //solve_stepwise_euler(f1,f2,1,prec);
  //solve_stepwise_heun(f1,f2,0.1,prec);
  check(limit(limit_test, REAL(1)), FUN1(1), prec);
  //solve_stepwise(f1,f2,0.1,100,prec);
  
  
  
}
