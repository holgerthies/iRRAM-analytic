#include "iRRAM.h"
#include "irram_analytic.h"
using namespace iRRAM;
using std::endl;

REAL sinseries1d(unsigned long n){
  if(n % 2 == 0) return 0;
  if((n-1) % 4 == 0) return inv_factorial(n);
  return -inv_factorial(n);
}
// sin(x)+1
REAL sinseriesmod1d(unsigned long n){
  if(n==0) return 1;
  if(n % 2 == 0) return 0;
  if((n-1) % 4 == 0) return inv_factorial(n);
  return -inv_factorial(n);
}

// series for sin(x1*x2*x3)
REAL sinseries(unsigned long n, unsigned long m, unsigned long q){
  if(n != m || n != q) return 0;
  if(n%2 == 0)
    return 0;
  else {
    if (0 == (n-1)%4)
      return inv_factorial(n);
    else
      return -inv_factorial(n);
  }
}

// series for sin(x1*x2*x3)
REAL sinseries2d(unsigned long n, unsigned long m){
  if(n != m ) return 0;
  if(n%2 == 0)
    return 0;
  else {
    if (0 == (n-1)%4)
      return inv_factorial(n);
    else
      return -inv_factorial(n);
  }
}
// series for sin(x1*x2*x3)+1
REAL sinseriesmod(unsigned long n, unsigned long m, unsigned long q){
  if(n+m+q==0) return 1;
  if(n != m || n != q) return 0;
  if(n%2 == 0)
    return 0;
  else {
    if (0 == (n-1)%4)
      return inv_factorial(n);
    else
      return -inv_factorial(n);
  }
}

struct sine : REAL_ANALYTIC<1>
{

  REAL get_coeff(const size_t n) const override
  {
    if(n % 2 == 0) return 0;
    if((n-1) % 4 == 0) return inv_factorial(n);
    return -inv_factorial(n);
  }

  REAL get_r() const override
  {
    return 2;
  }

  REAL get_M(const REAL& r) const override
  {
    return 6;
  }
      
};

struct sine3d : REAL_ANALYTIC<3>
{
  REAL get_coeff(const size_t n, const size_t m, const size_t q) const override
  {
    if(n != m || n != q) return 0;
    if(n%2 == 0)
      return 0;
    else {
      if (0 == (n-1)%4)
        return inv_factorial(n);
      else
        return -inv_factorial(n);
    }
  }

  REAL get_r() const override
  {
    return 2;
  }

  REAL get_M(const REAL& r) const override
  {
    return 6;
  }
      
};

void checkResult(const REAL& y,const REAL& sol, const int prec ){
  iRRAM::cout << "result: " << endl;
  rwrite(y, prec);
  iRRAM::cout << endl;
  iRRAM::cout << "should be " << endl;
  rwrite(sol,prec);
  iRRAM::cout << endl;
  if(!bound(abs(sol-y), -prec)){
    iRRAM::cout << "ERRROR" << endl;
  } else {
    iRRAM::cout << "OK!" << endl;
  }
}

void compute(){
  auto g = make_analytic<1,sine>();
  auto f = make_analytic<3,sine3d>();
  //auto sin_pwr = make_series<1,REAL>(std::function<REAL(const REAL&)>([] (const REAL& x) {return sin(x);}), 2,2);
	
  // auto sin_pwr3 = make_series<2,REAL>(std::function<REAL(const REAL&, const REAL&)>([] (const REAL& x, const REAL& y) {return sin(x*y);}), 2,2);
  //cout << get_coeff(*sin_pwr3, 1,1) << std::endl;
 //  auto f = make_analytic<REAL,REAL,REAL,REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(sinseries), 2, 2);
 //  auto g = make_analytic<REAL,REAL>(std::function<REAL(unsigned long)>(sinseries1d), 2,2);
 //  // ANALYTIC<1,REAL> h( 2,2, sin_pwr);
 //  auto h2 = make_analytic<REAL,REAL,REAL>(std::function<REAL(unsigned long, unsigned long)>(sinseries2d), 2, 2);
 //  //auto g = AnalyticFunction::projection<3,1,REAL>();
 //  //g = power(g,5);
  int l,prec;
  iRRAM::cin >>l>> prec;
  //  // rwrite(g({0.2,0.3, 0.5}), prec);
  //  iRRAM::cout << endl;
  //  // f continuation prec
  REAL x1= REAL(1)/REAL(4*l);
  REAL x2= REAL(1)/REAL(2*l);
  REAL x3= REAL(1)/REAL(8*l);

  REAL y = g->evaluate(x1);
  REAL sol=sin(x1);
  iRRAM::cout << "checking sin(x)" << endl;
  checkResult(y, sol, prec);
  
  auto x = variable_symbol<1,0>();
  auto h = sin(x);
  y = h->evaluate(x1);
  sol=sin(x1);
  iRRAM::cout << "checking predefined sin(x)" << endl;
  iRRAM::cout << h->to_string() << std::endl;
  checkResult(y, sol, prec);

  auto hd = pderive(h, 0, 102);
  iRRAM::cout << hd->to_string() << std::endl;
  simplify(hd);
  y = hd->evaluate(x1);
  sol=-sin(x1);
  iRRAM::cout << "checking predefined sin(x) derivative" << endl;
  iRRAM::cout << hd->to_string() << std::endl;
  checkResult(y, sol, prec);

  y = f->evaluate(x1,x2,x3);
  sol=sin(x1*x2*x3);
  iRRAM::cout << "checking sin(x1*x2*x3)" << endl;
  checkResult(y, sol, prec);

  auto sum = f+REAL(1)-f+f-REAL(2)+f-f+f;
  y = sum->evaluate(x1,x2,x3);
  iRRAM::cout << "checking addition" << endl;
  sol=2*sin(x1*x2*x3)-REAL(1);
  checkResult(y, sol, prec);

  auto prod1 = REAL(3)*g*g*REAL(2);
  y = prod1->evaluate(x1);
  iRRAM::cout << "checking multiplication 1d" << endl;
  sol=6*sin(x1)*sin(x1);
  checkResult(y, sol, prec);

  auto prod = REAL(3)*f*sum*REAL(2);
  y = prod->evaluate(x1,x2,x3);
  iRRAM::cout << "checking multiplication" << endl;
  sol=6*(2*sin(x1*x2*x3)-1)*sin(x1*x2*x3);
  checkResult(y, sol, prec);

  iRRAM::cout << "checking analytic continuation (1d)" << endl;
  // sin(x)/(sin(x)+1
  auto cont = transpose(g, {REAL(x1)});
  y = cont->evaluate(x1);
  sol=sin(x1+x1);
  checkResult(y, sol, prec);

  iRRAM::cout << "checking derivative (1d)" << endl;
  // sin(x)
  auto d1 = pderive(g*g,0,2);
  std::cout << d1->to_string() << std::endl;
  //d1 = d1->simplify();
  
  y = d1->evaluate(x1);
  //std::cout << d1->to_string() << std::endl;
  sol=2*cos(2*x1);
  
  checkResult(y, sol, prec);

 //  iRRAM::cout << "checking derivative (3d)" << endl;
 //  //auto d2 = f*pderive(f*pderive(f*f+f, 1, 1),2,100);
 //  auto d2 = pderive(pderive(pderive(f, 1, 1),2,1),0,1);
  
 //  std::cout << d2->to_string() << std::endl;
 //  //d2 = d2->simplify();
  
 // //  std::cout << d2->to_string() << std::endl;
 //  y = d2->evaluate(x1, x2,x3);
 //  //sol=-cos(x1*x2*x3)*x1*x1*x1*x2*x2*x3-sin(x1*x2*x3)*x1*x1*x2-x1*x1*x2*sin(x1*x2*x3);
 //  sol=2*x1*x2*x3*cos(2*x1*x2*x3)+sin(2*x1*x2*x3)+x1*(4*x2*x3*cos(2*x1*x2*x3)-4*x1*x2*x2*x3*x3*sin(2*x1*x2*x3));
 //  checkResult(y, sol, prec);

 //  // iRRAM::cout << "checking analytic continuation (3d)" << endl;
 //  // // sin(x)/(sin(x)+1
 //  // auto cont2 = continuation(f, 2,2, x3, x3, x3);
 //  // y = cont2->evaluate(x1,x2,x3);
 //  // sol=sin((x1+x3)*(x2+x3)*(x3+x3));
 //  // checkResult(y, sol, prec);

  auto comp = compose(g,g);
  iRRAM::cout << "checking composition (1d)" << endl;
  y = comp->evaluate(x1);
  sol=sin(sin(x1));
  checkResult(y, sol, prec);

  auto comp_deriv = pderive(comp, 0, 4);
  simplify(comp_deriv);
  //std::cout << comp_deriv->to_string() << std::endl;
  iRRAM::cout << "checking composition derivative (1d)" << endl;
  y = comp_deriv->evaluate(x1);
  sol=(1+6*cos(x1)*cos(x1))*cos(sin(x1))*sin(x1)+(-3+7*cos(x1)*cos(x1)+cos(x1)*cos(x1)*cos(x1)*cos(x1))*sin(sin(x1));
  checkResult(y, sol, prec);

  auto comp3 = compose(g,f);
  iRRAM::cout << "checking composition (3d)" << endl;
  y = comp3->evaluate(x1,x2,x3);
  sol=sin(sin(x1*x2*x3));
  checkResult(y, sol, prec);

  auto comp_deriv3 = pderive(pderive(comp3, 0, 1), 2,1);
  simplify(comp_deriv3);
  iRRAM::cout << "checking composition derivative (3d)" << endl;
  y = comp_deriv3->evaluate(x1,x2,x3);
  sol= x2*(cos(sin(x1*x2*x3))*(cos(x1*x2*x3)-x1*x2*x3*sin(x1*x2*x3))-x1*x2*x3*cos(x1*x2*x3)*cos(x1*x2*x3)*sin(sin(x1*x2*x3)));
  checkResult(y, sol, prec);
 //  auto comp4 = compose(g,h2)->to_analytic();
 //  iRRAM::cout << "checking composition (2d) by power series" << endl;
 //  y = comp4->evaluate(x1,x2);
 //  sol=sin(sin(x1*x2));
 //  checkResult(y, sol, prec);

 //  // auto cos1 = cos(h2);
 //  // iRRAM::cout << "checking cosine (2d)" << endl;
 //  // y = cos1->evaluate(x1,x2);
 //  // sol=cos(sin(x1*x2));
 //  // checkResult(y, sol, prec);

 //  // auto cos2 = cos(h2)->to_analytic();
 //  // iRRAM::cout << "checking cosine (2d) by power series" << endl;
 //  // y = cos2->evaluate(x1,x2);
 //  // sol=cos(sin(x1*x2));
 //  // checkResult(y, sol, prec);


 //  // auto exp1 = exp(h2);
 //  // iRRAM::cout << "checking exp (2d)" << endl;
 //  // y = exp1->evaluate(x1,x2);
 //  // sol=exp(sin(x1*x2));
 //  // checkResult(y, sol, prec);

 //  // auto exp2 = exp(h2)->to_analytic();
 //  // iRRAM::cout << "checking exp (2d) by power series" << endl;
 //  // y = exp2->evaluate(x1,x2);
 //  // sol=exp(sin(x1*x2));
 //  // checkResult(y, sol, prec);

  
 //  auto ssum = ((f+REAL(4))+REAL(5))+REAL(6)+f;
 //  std::cout << ssum->to_string() << std::endl;
  
 //  ssum = ssum->simplify();
 //  std::cout << ssum->to_string() << std::endl;
 //  y = ssum->evaluate(x1,x2,x3);
 //  iRRAM::cout << "checking scalar addition" << endl;
 //  sol=2*sin(x1*x2*x3)+REAL(15);
 //  checkResult(y, sol, prec);
  

 //  auto sum2 = (f-REAL(1)-f+f+REAL(2)+f)->to_analytic();
 //  y = sum2->evaluate(x1,x2,x3);
 //  iRRAM::cout << "checking addition by power series" << endl;
 //  sol=2*sin(x1*x2*x3)+REAL(1);
 //  checkResult(y, sol, prec);
   
 //  auto prod = REAL(3)*f*sum*REAL(2);
 //  y = prod->evaluate(x1,x2,x3);
 //  iRRAM::cout << "checking multiplication" << endl;
 //  sol=6*(2*sin(x1*x2*x3)+1)*sin(x1*x2*x3);
 //  checkResult(y, sol, prec);

 //  if(prec < 30)
 //  {
 //    auto prod2 = ( REAL(3)*f*sum*REAL(2) )->to_analytic();
 //    y = prod2->evaluate(x1,x2,x3);
 //    iRRAM::cout << "checking multiplication by power series" << endl;
 //    sol=6*(2*sin(x1*x2*x3)+1)*sin(x1*x2*x3);
 //    checkResult(y, sol, prec);
 //  }
 //  auto gmod = make_analytic<REAL,REAL>(std::function<REAL(unsigned long)>(sinseriesmod1d), 3,2);
 //  auto fmod = make_analytic<REAL,REAL,REAL,REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(sinseriesmod), 3,2);
  


 //  // iRRAM::cout << "checking division (1d)" << endl;
 //  // // sin(x)/(sin(x)+1)
 //  // auto inv_test = (g/gmod)/REAL( 2 );   
 //  // y = inv_test->evaluate(x1);
 //  // sol=sin(x1)/(1+sin(x1))/2;
 //  // checkResult(y, sol, prec);


 //  // iRRAM::cout << "checking division (1d) by power series" << endl;
 //  // // sin(x)/(sin(x)+1)
 //  // auto inv_test2 = ( (g/gmod)/REAL(2))->to_analytic();   
 //  // y = inv_test2->evaluate(x1);
 //  // sol=sin(x1)/(1+sin(x1))/2;
 //  // checkResult(y, sol, prec);


 //  // iRRAM::cout << "checking division (3d)" << endl;

 //  // // sin(x*y*z)/(sin(x*y*z)+1)
 //  // auto inv_test3d = f/fmod;   
 //  // y = inv_test3d->evaluate(x1,x2,x3);
 //  // sol=sin(x1*x2*x3)/(1+sin(x1*x2*x3));
 //  // checkResult(y, sol, prec);

 //  // iRRAM::cout << "checking division (3d) by power series" << endl;

 //  // // sin(x*y*z)/(sin(x*y*z)+1)
 //  // auto inv_test3d2 = (f/fmod)->to_analytic();   
 //  // y = inv_test3d2->evaluate(x1,x2,x3);
 //  // sol=sin(x1*x2*x3)/(1+sin(x1*x2*x3));
 //  // checkResult(y, sol, prec);

}
