#include "iRRAM.h"
#include "IVP.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;
using std::vector;

// 4y/t
REAL series(const unsigned long n, const unsigned long m){
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
	
  ANALYTIC<2,REAL> f(std::function<REAL(unsigned long, unsigned long)>(series), 0.25, 24);
  std::shared_ptr<IVP<2,REAL>> P(new IVP<2,REAL>({f}));


	int l,prec;
	iRRAM::cin >>l>> prec;
	// f continuation prec
	REAL x1= REAL(1)/REAL(16*l);
	REAL x2= REAL(1)/REAL(8*l);
  iRRAM::cout << "computing f("<<x1<<","<<x2<<")..."<<endl;
	REAL y = f(x1,x2);
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
  auto F = solve(P)[0];
  iRRAM::cout << "computing F("<<x1<<")..."<<endl;
	REAL Y = F(x1);
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
  iRRAM::cout << F.get_coeff(0)<<"+"<<F.get_coeff(1)<<"x+"<<F.get_coeff(2)<<"x^2" << endl;
  iRRAM::cout << F.get_coeff(3)<<"x^3+"<<F.get_coeff(4)<<"x^4+"<<F.get_coeff(5)<<"x^5" << endl;
}
