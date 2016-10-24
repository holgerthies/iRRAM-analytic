#include "iRRAM.h"
#include "irram_analytic.h"
#include "odes.h"
#include <chrono>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
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

template<class T>
REAL solve(const T& system, int method, DEBUG_INFORMATION& d )
{
  REAL x = 0.05;
  
  decltype(solve_taylor(system, x, d)) sols;
  if(method == 3)
    sols = solve_taylor(system,x,d);
  if(method == 4)
    sols = solve_taylor_deriv(system,x,d);
  return sols[1];
}

void compute(){
  DEBUG_INFORMATION d;
  
  vector<decltype(A2_SYSTEM())> systems_1d = {A2_SYSTEM(), A3_SYSTEM()};//, A5_SYSTEM(1.5)};
  vector<decltype(B1_SYSTEM())> systems_2d = {B1_SYSTEM()};
  
  int dimension, system,max_iter,method;
  struct rusage usage;
  struct timeval start, end;
  iRRAM::cin >> dimension >> system >> method >> max_iter;
  static int iteration_counter = 0;
  iteration_counter++; 
  if(iteration_counter == max_iter) return;
  
  
  getrusage(RUSAGE_SELF, &usage);
  start = usage.ru_utime;
  REAL sol;
  if(dimension == 1){
    sol = solve(systems_1d[system], method, d);
  }
  if(dimension == 2)
    sol = solve(systems_2d[system], method, d);
  sizetype error;
  sol.geterror(error);
  getrusage(RUSAGE_SELF, &usage);
  end = usage.ru_utime;
  auto iteration_time =  end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1000000;
  int error_exp_normalized;
  unsigned long mantissa = error.mantissa;
  error_exp_normalized = error.exponent;
  while(mantissa > 1){
    mantissa /= 2;
    error_exp_normalized++;
  }
  
  check(sol, A3_sol(0.05), 100);
  return;
  
  
  std::cout << " dimension:" << dimension;
  std::cout << " system:" << system;
  std::cout << " method:" << method;
  std::cout << " time:" << iteration_time;
  std::cout << " steps:" << d.steps;
  std::cout << " order:" << d.order;
  std::cout << " iteration:" << iteration_counter;
  std::cout << " precision:" << ACTUAL_STACK.actual_prec;
  std::cout << " error_mantissa:"  << error.mantissa;
  std::cout << " error_exponent:" << error.exponent;
  std::cout << " normalized:" << error_exp_normalized;
  
  
  std::cout << std::endl;
  
  
  iRRAM::cout << (REAL(0) == REAL(0)) << endl;
  
  
  
}
