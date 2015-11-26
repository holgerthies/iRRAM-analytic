#include "iRRAM.h"
#include "POWERSERIES.h"
#include <iostream>
#include <stdlib.h> 
#include <time.h>
using namespace std;
using namespace iRRAM;

template<unsigned int n, typename... ARGS>
REAL check_poly(const POLYNOMIAL<n,REAL>& p, const REAL& x, const ARGS&... rest){
  REAL ans=0;
  REAL pow=1;
  for(int i=0; i<p.get_degree(); i++){
    ans += pow*check_poly(p.get_coefficient(i), rest...); 
    pow *= x;
  }
  return ans;
}

template<>
REAL check_poly(const POLYNOMIAL<1,REAL>& p, const REAL& x){
  REAL ans=0;
  REAL pow=1;
  for(int i=0; i<p.get_degree(); i++){
    ans += pow*p.get_coefficient(i); 
    pow *= x;
  }
  return ans;
}
POLYNOMIAL<1,REAL> rand_poly1d(){
  srand (time(NULL));
  int deg=rand() % 1000;
  std::vector<REAL> c(deg);
  for(int i=0; i<deg;i++){
    c[i] = REAL(rand() % 1000)/REAL(rand() % 1000 + 1);
  }
  return  POLYNOMIAL<1,REAL>(c);
}

template<unsigned int n>
POLYNOMIAL<n,REAL> rand_poly(){
  int deg=rand() % 10;
  std::vector<POLYNOMIAL<n-1,REAL>> c;
  for(int i=0; i<deg;i++){
    c.push_back(rand_poly<n-1>());
  }
  return  POLYNOMIAL<n,REAL>(c);
}

template<>
POLYNOMIAL<1,REAL> rand_poly<1>(){
  return rand_poly1d();
}

void compute(){
  srand (time(NULL));
  REAL x = REAL(rand() % 1000)/REAL(rand() % 1000 + 1);
  REAL y = REAL(rand() % 1000)/REAL(rand() % 1000 + 1);
  iRRAM::cout << x <<","<<y<<std::endl;
  iRRAM::cout << std::endl;
  auto test = rand_poly1d();
  auto test2 = rand_poly1d();
  auto test2d = rand_poly<2>();
  auto test7d = rand_poly<7>();
  
  iRRAM::cout << test(x) << std::endl;
  iRRAM::cout << check_poly(test,x) << std::endl;
  iRRAM::cout << std::endl;
  iRRAM::cout << test2(x) << std::endl;
  iRRAM::cout << check_poly(test2,x) << std::endl;
  iRRAM::cout << std::endl;
  std::cout << get_total_degree(test2d)<<std::endl;
  iRRAM::cout << test2d(x,y) << std::endl;
  iRRAM::cout << check_poly(test2d, x,y) << endl;
  std::cout << get_total_degree(test7d)<<std::endl;
  iRRAM::cout << test7d(x,x,x,y,y,x,y) << std::endl;
  iRRAM::cout << check_poly(test7d,x,x,x,y,y,x,y) << endl;
}
