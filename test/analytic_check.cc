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
REAL sinseries(const std::vector<unsigned long>& v){
  int n=v[0];
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
	
  ANALYTIC f(std::shared_ptr<std::function<REAL(const std::vector<unsigned long>&)>>(new std::function<REAL(const std::vector<unsigned long>&)>(sinseries)), 2);
	int l,prec;
	iRRAM::cin >>l>> prec;
	// f continuation prec
	REAL x= REAL(1)/REAL(4*l);
  REAL x_real = x;
	REAL y = f.eval_k(x, 0); 
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
