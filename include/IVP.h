#ifndef _IVP_H
#define _IVP_H 
#include "ANALYTIC.h"
#include <vector>
namespace iRRAM{
  template <unsigned int d, class T>
  class IVP{
      private:
        std::shared_ptr<ANALYTIC<d,T>> f;
        // precomputed values for the coefficients in the recursion formula
        std::vector<REAL> a;
        // ai[i,n] contains the coefficient a(i,i+n) since a(i,n)=0 for i>n
        std::vector<std::vector<REAL>> ai;
        REAL get_a(const int n, const int i);
      public:
        IVP(const std::shared_ptr<ANALYTIC<d,T>> f ): f(f) {};
        IVP(const ANALYTIC<d,T> f): IVP(std::shared_ptr<ANALYTIC<d,T>>(new ANALYTIC<d,T>(f))) {};
        REAL get_coeff(const int n);
        REAL get_r() const {return f->get_r();}
        REAL get_M() const {return f->get_r()*f->get_M();}
  };
  
  template<unsigned int d, class T>
  REAL IVP<d,T>::get_a(const int i, const int n){
    if(i > n) return 0;
    int m=n-i;
    if(ai.size() > i && ai[i].size() > m) return ai[i][m];
    if(ai.size() <= i){
      if(i>0)
        get_a(i-1, n); // make sure ai.size() = i
      // add empty vector at position i
      ai.push_back(std::vector<REAL>());
    }
    if(ai[i].size() <= m){
      if(m > 0)
        get_a(i, m-1); // make sure ai[i].size() = m
      ai[i].push_back(0);
    }
    if(i == 0 && n==0) ai[i][m] = 1;
    else if(i > 0){
      for(unsigned long j=0; j<=n;j++){
          ai[i][m]+=get_a(i-1,n-j)*get_coeff(j);
      } 
    }
    return ai[i][m];
  }

  template<unsigned int d, class T>
  REAL IVP<d,T>::get_coeff(const int n){
    if(a.size() > n) return a[n];
    // make sure a.size() = n
    if(n > 0)
      get_coeff(n-1);
    a.push_back(0);
    if(n > 0)
    {
      for(auto w : partitions(n-1, 2)){
        auto l = w[0];
        auto k = w[1];
        for(unsigned long i=0; i<=l; i++)
        {
          a[n] += f->get_coeff({k,i})*get_a(i,l);
        }
      }
      a[n] /= REAL(n);
    }
    return a[n];
  };

  template <class T>
  ANALYTIC<1,T> solve(std::shared_ptr<IVP<2,T>> P){
    auto series = [P] (const std::vector<unsigned long>& v) {
      return P->get_coeff((int)v[0]);
    };
    return ANALYTIC<1,T>(series, P->get_r(), P->get_M());
  }
}
#endif
