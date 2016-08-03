#ifndef ODES_H
#define ODES_H
#include "iRRAM.h"
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

REAL A2_sol(const REAL& t)
{
  REAL C=1;
  return REAL(1)/sqrt(t+C);
}


// Problem A3
// an oscillatory problem
REAL A3_fun(const REAL& t, const REAL& y)
{
  return y*cos(t);
}

REAL A3_sol(const REAL& t)
{
  REAL C=1;
  return C*exp(sin(x));
}

// Problem A5
// a spiral curve
REAL A5_fun(const REAL& t, const REAL& y)
{
  return (y-t)/(y+t);
}

// Problem Class B: Small systems

// Problem B1
// the growth of two confilicting populations
REAL B1_fun1(const REAL& t, const REAL& y1, const REAL& y2)
{
  return 2*(y1-y1*y2);
}
REAL B1_fun2(const REAL& t, const REAL& y1, const REAL& y2)
{
  return -(y2-y1*y2);
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
return 0.01-(1+(y1+1000)*(y1+1))*(0.01+y1+y2);
}
REAL SD5_fun2(const REAL& t, const REAL& y1, const REAL& y2)
{
return 0.01-(1+y2*y2)*(0.01+y1+y2);
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

