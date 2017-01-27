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
REAL solve(const T& system, const REAL& x, const int method, const int solver, const bool out, DEBUG_INFORMATION& d )
{
  
  decltype(ivp_solve_cs(system, x,out,0, d)) sols;
  if(method == 1)
    sols = ivp_solve_cs(system,x,out,solver, d);
  if(method == 2)
    sols = ivp_solve_co(system,x,25, out, d);
  if(method == 3)
    sols = ivp_solve_co(system,x,50, out, d);
  if(method == 4)
    sols = ivp_solve_mixed(system,x, out, d);
  return sols[0];
}


void compute(){
  DEBUG_INFORMATION d;
  static int iteration_counter = 0;

  vector<decltype(A2_SYSTEM())> systems_1d = {A2_SYSTEM()};//, A3_SYSTEM(), A5_SYSTEM(), SINCOS_SYSTEM()};
  vector<decltype(A3_SYSTEM())> systems_2d = {A3_SYSTEM(), E2_SYSTEM()};
  vector<decltype(SINCOS_POLY_SYSTEM())> systems_6d = {SINCOS_POLY_SYSTEM()};
  
  int dimension=0, system,max_iter,prec,method,solver_type;
  struct rusage usage;
  struct timeval start, end;
  
  while(dimension == 0){
    iRRAM::cout << "choose system dimension" << std::endl;
    iRRAM::cin >> dimension;
    if(dimension != 1 && dimension != 2 && dimension != 6){
      iRRAM::cout << "invalid dimension" << std::endl;
      dimension = 0;
    }
  }
  iRRAM::cout << "choose system" << std::endl;
  iRRAM::cin >>  system;
  
  iRRAM::cout << "choose solver method" << std::endl;
  iRRAM::cin >>  solver_type;

  iRRAM::cout << "choose step size method" << std::endl;
  iRRAM::cin >>  method;

  REAL x;
  iRRAM::cout << "choose x" << std::endl;
  iRRAM::cin >>  x;

  iRRAM::cout << "choose precision (or 0 for iteration number)" << std::endl;
  iRRAM::cin >>  prec;
  bool out = true;
  if(prec > 0){
    iRRAM::cout << setRwidth(prec) << std::endl;
  }
  else{
    iRRAM::cout << "choose number of iterations" << std::endl;
    iRRAM::cin >>  max_iter;
    out=false;
  }

  iteration_counter++; 
  if(prec <= 0 && iteration_counter == max_iter) return;
  
  
  getrusage(RUSAGE_SELF, &usage);
  start = usage.ru_utime;
  REAL sol;
  
  if(dimension == 1){
    sol = solve(systems_1d[system],x, method,solver_type, out,  d);
  }
  if(dimension == 2)
    sol = solve(systems_2d[system],x, method,solver_type, out,  d);
  if(dimension == 6)
    sol = solve(systems_6d[system],x, method,solver_type, out,  d);
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
  
  //check(sol, A3_sol(5), 50);
  // return;
  
  
  std::cout << " dimension:" << dimension;
  std::cout << " system:" << system;
  std::cout << " time:" << iteration_time;
  std::cout << " steps:" << d.steps;
  std::cout << " order:" << d.order;
  std::cout << " iteration:" << iteration_counter;
  std::cout << " precision:" << ACTUAL_STACK.actual_prec;
  std::cout << " error_mantissa:"  << error.mantissa;
  std::cout << " error_exponent:" << error.exponent;
  std::cout << " normalized:" << error_exp_normalized;
  
  
  std::cout << std::endl;
  
  if(prec <= 0)
    iRRAM::cout << (REAL(0) == REAL(0)) << endl;
  
  
  
}
