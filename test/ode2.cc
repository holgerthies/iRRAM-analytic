#include "iRRAM.h"
#include "IVP.h"
#include "combinatorics.h"
using namespace iRRAM;
using std::endl;
using std::vector;
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

// f(x,y,z) = z
REAL series1(const std::vector<unsigned long>& v){
  if(v[0] == 0 && v[1] == 0 && v[2] == 0) return 1;
  if(v[0] == 0 && v[1] == 0 && v[2] == 1) return 1;
  return 0;
}

// f(x,y,z) = y
REAL series2(const std::vector<unsigned long>& v){
  if(v[0] == 0 && v[1] == 1 && v[2] == 0) return -1;
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
void compute(){
	
  ANALYTIC<3,REAL> f1(std::shared_ptr<std::function<REAL(const std::vector<unsigned long>&)>>(new std::function<REAL(const std::vector<unsigned long>&)>(series1)),1, 2);
  ANALYTIC<3,REAL> f2(std::shared_ptr<std::function<REAL(const std::vector<unsigned long>&)>>(new std::function<REAL(const std::vector<unsigned long>&)>(series2)),1, 2);
  std::shared_ptr<IVP<3,REAL>> P(new IVP<3,REAL>({f1,f2}));

	int l,prec;
	iRRAM::cin >>l>> prec;
	// f continuation prec
	REAL x1= REAL(1)/REAL(16*l);
	REAL x2= REAL(1)/REAL(8*l);
	REAL x3= REAL(1)/REAL(10*l);
  iRRAM::cout << "computing f1("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
	REAL y = f1({x1,x2,x3});
	iRRAM::cout << "result: " << endl;
	rwrite(y, prec);
	iRRAM::cout << endl;
	REAL sol=fun1(x1,x2,x3);
	iRRAM::cout << "should be " << endl;
	rwrite(sol,prec);
	iRRAM::cout << endl;
	if(!bound(abs(sol-y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
  iRRAM::cout << "computing f2("<<x1<<","<<x2<<","<< x3 << ")..."<<endl;
	y = f2({x1,x2,x3});
	iRRAM::cout << "result: " << endl;
	rwrite(y, prec);
	iRRAM::cout << endl;
	sol=fun2(x1,x2,x3);
	iRRAM::cout << "should be " << endl;
	rwrite(sol,prec);
	iRRAM::cout << endl;
	if(!bound(abs(sol-y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
  iRRAM::cout << "solving ode..." << endl;
  auto F = solve(P);
  iRRAM::cout << "computing F1("<<x1<<")..."<<endl;
	REAL Y = F[0]({x1});
	iRRAM::cout << "result: " << endl;
	rwrite(Y, prec);
	iRRAM::cout << endl;
	REAL SOL=FUN1(x1);
	iRRAM::cout << "should be " << endl;
	rwrite(SOL,prec);
	iRRAM::cout << endl;
	if(!bound(abs(SOL-Y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
  iRRAM::cout << "computing F2("<<x1<<")..."<<endl;
	Y = F[1]({x1});
	iRRAM::cout << "result: " << endl;
	rwrite(Y, prec);
	iRRAM::cout << endl;
	SOL=FUN2(x1);
	iRRAM::cout << "should be " << endl;
	rwrite(SOL,prec);
	iRRAM::cout << endl;
	if(!bound(abs(SOL-Y), -prec)){
		iRRAM::cout << "ERRROR" << endl;
	} else {
		iRRAM::cout << "OK!" << endl;
	}
}