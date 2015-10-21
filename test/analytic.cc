#include "iRRAM.h"
#include "ANALYTIC.h"
#include <ctime>
#include <cstdlib>
#include <map>
#include <chrono>
#include <ctime>
#include <sys/time.h>
#include <sys/resource.h>
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
REAL sinseries(const int n){
	if (0 == n%2)
		return 0;
	else {
		if (0 == (n-1)%4)
			return inv_factorial(n);
		else
			return -inv_factorial(n);
	 }
}

void compute(){
	
	int l, max_iter, continuation_num;
  iRRAM::cin >> l >> continuation_num >> max_iter;
	static int iteration_counter = 0;
	iteration_counter++; 
  if(iteration_counter == max_iter) return;
	// f continuation prec
	REAL x = REAL(1)/REAL(4*l);
  struct rusage usage;
  struct timeval start, end;
  ANALYTIC f(std::shared_ptr<std::function<REAL(unsigned long)>>(new std::function<REAL(unsigned long)>(sinseries)), 2);
	sizetype error;
  getrusage(RUSAGE_SELF, &usage);
  start = usage.ru_utime;
	REAL y = f.eval_k(x,continuation_num); 
	y.geterror(error);
  getrusage(RUSAGE_SELF, &usage);
  end = usage.ru_utime;
  auto iteration_time =  end.tv_sec-start.tv_sec+double(end.tv_usec-start.tv_usec)/1000000;
  std::cout << "l:" << l;
  std::cout << " continuations:" << continuation_num;
  std::cout << " time:" << iteration_time;
	std::cout << " iteration:" << iteration_counter;
	std::cout << " precision:" << ACTUAL_STACK.actual_prec;
	std::cout << " error_mantissa:"  << error.mantissa;
	std::cout << " error_exponent:" << error.exponent;
	int numseries = max(1, 2*continuation_num);
	std::cout << " cached_coefficients: ";
	for(int i=0; i<numseries;i++){
		std::cout << f.get_known_coeffs(i) << " ";
	}
	std::cout << endl;
	//iRRAM::cout << "result" << endl;
  iRRAM::cout << (REAL(0) == REAL(0)) << endl;
	iRRAM::cout << "end" << endl;
	
}
