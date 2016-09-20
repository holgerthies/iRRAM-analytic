#ifndef ODES_H
#define ODES_H

#include "iRRAM.h"
#include "ANALYTIC.h"
#include "IVPSOLVER.h"
#include "combinatorics.h"
using namespace iRRAM;

  
// Class A: Single equations
// Problem A2
// A special case of the Riccati Equation
REAL A2_fun(const REAL& t, const REAL& y)
{
  return -power(y,3)/2;
}
REAL A2_series(const unsigned long v0, const unsigned long v1)
{
  if(v1 == 3 && v0==0) return -REAL(1)/REAL(2);
  return 0;
}

REAL A2_max(const REAL& r)
{
  
  return power(r,3)/2;
}

IVPSYSTEM<REAL, REAL, REAL> A2_SYSTEM(const REAL& r)
{

  auto A2_analytic = real_analytic<2>::make(std::function<REAL(unsigned long, unsigned long)>(A2_series),A2_max(r),r);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL>>(A2_analytic)->add_algorithm(A2_fun);
  IVPSYSTEM<REAL, REAL, REAL> A2_IVP;
  A2_IVP.F = {A2_analytic};
  A2_IVP.y = {0, 1};
  return A2_IVP;


}

REAL A2_sol(const REAL& t)
{
  REAL C=1;
  return REAL(1)/sqrt(t+C);
}


// // Problem A3
// // an oscillatory problem
REAL A3_fun(const REAL& t, const REAL& y)
{
  return y*cos(t);
}

REAL A3_series(const unsigned long v0, const unsigned long v1)
{
  
  if(v1 != 1) return 0;
  if(v0 % 2 == 1) return 0;
  if (v0 % 4 == 0)
      return inv_factorial(v0);
  return -inv_factorial(v0);
}

REAL A3_max(const REAL& r)
{
  return r*exp(r);
}

IVPSYSTEM<REAL, REAL, REAL> A3_SYSTEM(const REAL& r)
{

  auto A3_analytic = make_analytic<REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long)>(A3_series),A3_max(r),r);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL>>(A3_analytic)->add_algorithm(A3_fun);
  IVPSYSTEM<REAL, REAL, REAL> A3_IVP;
  A3_IVP.F = {A3_analytic};
  A3_IVP.y = {0, 1};
  return A3_IVP;


}

REAL A3_sol(const REAL& t)
{
  REAL C=1;
  return C*exp(sin(t));
}

// Problem A5
// a spiral curve
REAL A5_fun(const REAL& t, const REAL& y)
{
  return (y-t+4)/(y+t+4);
}

REAL A5_series(const unsigned long v0, const unsigned long v1)
{
  if(v0 == 0 && v1==0) return 0;
  REAL ans=8*REAL(factorial(v0+v1-1))*inv_factorial(v0-1)*inv_factorial(v1);
  if((v0+v1) % 2 == 0) return ans;
  return -ans;
}

REAL A5_max(const REAL& r)
{
  return (2*r+4)/(-2*r+4);
}

IVPSYSTEM<REAL, REAL, REAL> A5_SYSTEM(const REAL& r)
{

  auto A5_analytic = make_analytic<REAL,REAL, REAL>(std::function<REAL(unsigned long, unsigned long)>(A5_series),A5_max(r),r);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL>>(A5_analytic)->add_algorithm(A5_fun);
  IVPSYSTEM<REAL, REAL, REAL> A5_IVP;
  A5_IVP.F = {A5_analytic};
  A5_IVP.y = {0, 0};
  return A5_IVP;
}

// Problem Class B: Small systems

// Problem B1
// the growth of two confilicting populations
REAL B1_fun1(const REAL& t, const REAL& y1, const REAL& y2)
{
  return 2*(y1-y1*y2);
}

REAL B1_series1(const unsigned long v0, const unsigned long v1, const unsigned long v2)
{
  if(v1 == 1 && v2==0) return 2;
  if(v1 == 1 && v2==1) return -2;
  return 0;
}

REAL B1_max1(const REAL& r)
{
  return 2*r+r*r;
}

REAL B1_fun2(const REAL& t, const REAL& y1, const REAL& y2)
{
  return -(y2-y1*y2);
}

REAL B1_max2(const REAL& r)
{
  return r+r*r;
}

REAL B1_series2(const unsigned long v0, const unsigned long v1, const unsigned long v2)
{
  if(v1 == 0 && v2==1) return -1;
  if(v1 == 1 && v2==1) return 1;
  return 0;
}

IVPSYSTEM<REAL, REAL, REAL, REAL> B1_SYSTEM(const REAL& r)
{

  auto B1_analytic1 = make_analytic<REAL,REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(B1_series1),B1_max1(r),r);
  auto B1_analytic2 = make_analytic<REAL,REAL, REAL,REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(B1_series2),B1_max2(r),r);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(B1_analytic1)->add_algorithm(B1_fun1);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(B1_analytic2)->add_algorithm(B1_fun2);
  IVPSYSTEM<REAL, REAL, REAL, REAL> B1_IVP;
  B1_IVP.F = {B1_analytic1, B1_analytic2};
  B1_IVP.y = {0, 1, 3};
  return B1_IVP;
}

// Problem B3
// a nonlinear chemical reaction
REAL B3_fun1(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
  return -y1;
}
REAL B3_fun2(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
  return y1-y2*y2;
}
REAL B3_fun3(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
  return y2*y2;
}

// Problem B4
// the integral surface of a torus
REAL B4_fun1(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
  return -y2-y1*y3/(sqrt(y1*y1+y2*y2));
}
REAL B4_fun2(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
  return y1-y2*y3/(sqrt(y1*y1+y2*y2));
}
REAL B4_fun3(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
  return y1/(sqrt(y1*y1+y2*y2));
}

// Problem Class C: Moderate Problems

// Problem C1

template<int I, class... Args>
std::function<REAL(const Args&...)> C1_fun()
{
  using std::get;
  return [] (const Args&... args) -> REAL {
    auto args_t = std::make_tuple<const Args&...>(args...);
    if(I == 1) return -get<I>(args_t);
    if(I == sizeof...(Args)-1) return get<I-1>(args_t);
    return get<I-1>(args_t)-get<I>(args_t);
  };
}

unsigned long sum()
{
  return 0;
}

template<class... Args>
unsigned long sum(const unsigned long x, const Args... args)
{
  return x+sum(args...);
}

template<int I, class... Args>
std::function<REAL(Args...)> C1_series()
{
  using std::get;
  return [] (Args... args) -> REAL {
    auto args_t = std::make_tuple<Args...>(std::forward<Args>(args)...);
    
    if(sum(args...) != 1)
        return 0;
    if(I != 1 && get<I-1>(args_t) == 1 && get<I>(args_t) == 0)
      return 1;
    if(I != sizeof...(args)-1 && get<I-1>(args_t) == 0 && get<I>(args_t) == 1)
      return -1;
    return 0;
  };
  
  
}

template<int I, int D, class... Args>
struct C1_fun_getter
{
  using type = typename C1_fun_getter<I, D-1, Args..., REAL>::type;
  static type get()
  {
    return C1_fun_getter<I, D-1, Args..., REAL>::get();
  }
};


template<int I, class... Args>
struct C1_fun_getter<I,0, Args...>
{
  using type = decltype(C1_fun<I, Args...>());
  static type get()
  {
    return C1_fun<I, Args...>();
  }
};
template<int I, int D, class... Args>
struct C1_series_getter
{
  using type = typename C1_series_getter<I, D-1, Args..., unsigned long>::type;
  static type get()
  {
    return C1_series_getter<I, D-1, Args..., unsigned long>::get();
  }
};


template<int I, class... Args>
struct C1_series_getter<I,0, Args...>
{
  using type = decltype(C1_series<I, Args...>());
  static type get()
  {
    return C1_series<I, Args...>();
  }
};


REAL C1_max(const REAL& r)
{
  return 2*r;
}

template<int I, class... Args>
struct C1_analytic_getter
{
  using analytic_type = typename C1_analytic_getter<I-1, REAL, Args...>::analytic_type;
  using vector_type = typename C1_analytic_getter<I-1, REAL, Args...>::vector_type;
  using system_type = typename C1_analytic_getter<I-1, REAL, Args...>::system_type;
  
  static vector_type get(const REAL& r)
  {
    auto F = real_analytic<I+sizeof...(Args)+1>::make(C1_series_getter<I, I+sizeof...(Args)+1>::get(),C1_max(r),r);
    std::static_pointer_cast<analytic_type>(F)->add_algorithm(C1_fun_getter<I, I+sizeof...(Args)+1>::get());
    auto C1_analytic = C1_analytic_getter<I-1, REAL, Args...>::get(r);
    C1_analytic.push_back(F);
    return C1_analytic;
  }

};

template<class... Args>
struct C1_analytic_getter<0, Args...>
{
  using analytic_type = ANALYTIC<REAL, REAL, Args...>;
  using vector_type = std::vector<std::shared_ptr<Node<REAL, REAL, Args...>>>;
  using system_type = IVPSYSTEM<REAL, REAL, Args...>;
  
  static vector_type get(const REAL& r)
  {
    return vector_type();
  }
};


  
template<int dim>
typename C1_analytic_getter<dim-1>::system_type C1_SYSTEM(const REAL& r)
{

  auto C1_analytic = C1_analytic_getter<dim-1>::get(r);
  typename C1_analytic_getter<dim-1>::system_type C1_IVP;
  C1_IVP.F = C1_analytic;
  C1_IVP.y = std::vector<REAL>(dim, 0);
  C1_IVP.y[1] = 1;
  
  return C1_IVP;
}
//Problem Class D: Orbit Equations

//Problem D1
REAL D1_fun1(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -y3;
}
REAL D1_fun2(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return y4;
}
REAL D1_fun3(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -y1/power(y1*y1+y2*y2, REAL(3)/REAL(2));
}
REAL D1_fun4(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -y2/power(y1*y1+y2*y2, REAL(3)/REAL(2));
}

//Problem Class E: Higher Order Equations

// Problem E3
// Derived from Duffing's equation
REAL E3_fun1(const REAL& t, const REAL& y1, const REAL& y2)
{
  return y2;
}
REAL E3_fun2(const REAL& t, const REAL& y1, const REAL& y2)
{
  return power(y1,3)/REAL(6)-y1+2*sin(2.78535*t);
}


// Stiff Problems

// Problem Class A: Linear with Real Eigenvalues

// Problem A1
REAL SA1_fun1(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -0.5*y1;
}

REAL SA1_fun2(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -y2;
}
  
REAL SA1_fun3(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -100*y3;
}

REAL SA1_fun4(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
  return -90*y4;
}

// Problem Class B: Linear with non-real eigenvalues

// Problem B1
REAL SB1_fun1(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
return -y1+y2;
}

REAL SB1_fun2(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
return -100*y1-y2;
}
  
REAL SB1_fun3(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
return -100*y3+y4;
}

REAL SB1_fun4(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3, const REAL& y4)
{
return -10000*y3-100*y4;
}

// Problem Class D: Non-linear with real eigenvalues

// Problem D1
// nuclear reactor theory
REAL SD1_fun1(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
return REAL(0.2)*(y2-y1);
}
REAL SD1_fun2(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
return 10*y1-(60-0.125*y3)*y2+0.125*y3;
}
REAL SD1_fun3(const REAL& t, const REAL& y1, const REAL& y2, const REAL& y3)
{
return 1;
}
// Problem D5
// reactor kinetics
REAL SD5_fun1(const REAL& t, const REAL& y1, const REAL& y2)
{
  return REAL(0.01)-(1+(y1+1000)*(y1+1))*(REAL(0.01)+y1+y2);
}
REAL SD5_fun2(const REAL& t, const REAL& y1, const REAL& y2)
{
  return REAL(0.01)-(1+y2*y2)*(REAL(0.01)+y1+y2);
}

REAL SD5_max1(const REAL& r)
{
  return REAL(0.01)+(1+(r+1000)*(r+1))*(REAL(0.01)+2*r);
}

REAL SD5_max2(const REAL& r)
{
  return REAL(0.01)+(1+r*r)*(REAL(0.01)+2*r);
}

REAL SD5_series1(const unsigned long v0, const unsigned long v1, const unsigned long v2)
{
  if(v1 == 0 && v2==0) return -10;
  if(v1 == 1 && v2==0) return -1011.01;
  if(v1 == 0 && v2==1) return -1001;
  if(v1 == 1 && v2==1) return -1001;
  if(v1 == 2 && v2==0) return -1001.10;
  if(v1 == 2 && v2==1) return -1;
  if(v1 == 3 && v2==0) return -1;
  return 0;
}

REAL SD5_series2(const unsigned long v0, const unsigned long v1, const unsigned long v2)
{
  if(v1 == 1 && v2==0) return -1;
  if(v1 == 0 && v2==1) return -1;
  if(v1 == 0 && v2==2) return -0.01;
  if(v1 == 1 && v2==2) return -1;
  if(v1 == 0 && v2==3) return -1;
  return 0;
}

IVPSYSTEM<REAL, REAL, REAL, REAL> SD5_SYSTEM(const REAL& r)
{

  auto SD5_analytic1 = make_analytic<REAL,REAL, REAL, REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(SD5_series1),SD5_max1(r),r);
  auto SD5_analytic2 = make_analytic<REAL,REAL, REAL,REAL>(std::function<REAL(unsigned long, unsigned long, unsigned long)>(SD5_series2),SD5_max2(r),r);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(SD5_analytic1)->add_algorithm(SD5_fun1);
  std::static_pointer_cast<ANALYTIC<REAL, REAL, REAL, REAL>>(SD5_analytic2)->add_algorithm(SD5_fun2);
  IVPSYSTEM<REAL, REAL, REAL, REAL> SD5_IVP;
  SD5_IVP.F = {SD5_analytic1, SD5_analytic2};
  SD5_IVP.y = {0, 1, 3};
  return SD5_IVP;
}
// Problem Class E: non-linear with non real eigenvalues

// Problem E2
// Van der Pol's equation
REAL SE2_fun1(const REAL& t, const REAL& y1, const REAL& y2)
{
return y2;
}
REAL SE2_fun2(const REAL& t, const REAL& y1, const REAL& y2)
{
return 5*(1-y1*y1)*y2-y1;
}


#endif /* ODES_H */

