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
	
  ANALYTIC f(std::shared_ptr<std::function<REAL(unsigned long)>>(new std::function<REAL(unsigned long)>(sinseries)), 2);
	int l,prec, continuation_num;
	iRRAM::cin >>l>> prec >> continuation_num;
	// f continuation prec
	REAL x= REAL(1)/REAL(4*l);
  REAL x_real = x+REAL(continuation_num)/REAL(2*l);
	REAL y = f.eval_k(x, continuation_num); 
	iRRAM::cout << "result: " << endl;
	rwrite(y, prec);
	iRRAM::cout << endl;
	REAL sol=sin(x_real);
	iRRAM::cout << "should be " << endl;
	rwrite(sol,prec);
	iRRAM::cout << endl;
	if(!bound(abs(sol-y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
	
}
