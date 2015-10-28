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

// 4y/t
REAL series(const std::vector<unsigned long>& v){
  int n=v[0];
  int m=v[1];
  if(m >= 2) return 0;
  REAL ans=4;
  if(n % 2 == 0) return ans;
  return -ans;
}

REAL fun(const REAL& t, const REAL& y){
  return 4*(y+1)/(t+1);
}

REAL FUN(const REAL& t){
  return power(t+1, 4)-1;
}

void compute(){
	
  ANALYTIC<2,REAL> f(std::shared_ptr<std::function<REAL(const std::vector<unsigned long>&)>>(new std::function<REAL(const std::vector<unsigned long>&)>(series)), 0.25, 24);
	int l,prec;
	iRRAM::cin >>l>> prec;
	// f continuation prec
	REAL x1= REAL(1)/REAL(16*l);
	REAL x2= REAL(1)/REAL(8*l);
  iRRAM::cout << "computing f("<<x1<<","<<x2<<")..."<<endl;
	REAL y = f({x1,x2});
	iRRAM::cout << "result: " << endl;
	rwrite(y, prec);
	iRRAM::cout << endl;
	REAL sol=fun(x1,x2);
	iRRAM::cout << "should be " << endl;
	rwrite(sol,prec);
	iRRAM::cout << endl;
	if(!bound(abs(sol-y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
  iRRAM::cout << "solving ode..." << endl;
  auto F = solve(f);
  iRRAM::cout << "computing F("<<x1<<")..."<<endl;
	REAL Y = F({x1});
	iRRAM::cout << "result: " << endl;
	rwrite(Y, prec);
	iRRAM::cout << endl;
	REAL SOL=FUN(x1);
	iRRAM::cout << "should be " << endl;
	rwrite(SOL,prec);
	iRRAM::cout << endl;
	if(!bound(abs(SOL-Y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
  iRRAM::cout << "first coefficients " << endl;
  iRRAM::cout << F.get_coeff({0})<<"+"<<F.get_coeff({1})<<"x+"<<F.get_coeff({2})<<"x^2" << endl;
  iRRAM::cout << F.get_coeff({3})<<"x^3+"<<F.get_coeff({4})<<"x^4+"<<F.get_coeff({5})<<"x^5" << endl;
}