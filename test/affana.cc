#include "iRRAM.h"
#include "ANALYTIC.h"
#include "ADDITION.h"
#include "SUBTRACTION.h"
#include "MULTIPLICATION.h"
#include "DIVISION.h"
#include "DERIVATIVE.h"
#include "COMPOSITION.h"
#include "IVPSOLVER.h"
#include "CONTINUATION.h"
#include "SINE.h"
#include "EXPONENTIATION.h"
#include "coefficient_computation.h"
#include "combinatorics.h"
#include <iostream>
using namespace std;
using namespace iRRAM;
// typedef REAL RTYPE;
REAL sinseries1d(unsigned long n){
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

void compute(){
  auto g = make_analytic<REAL,REAL>(std::function<REAL(unsigned long)>(sinseries1d), 2,2);
}

