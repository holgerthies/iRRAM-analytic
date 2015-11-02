#ifndef _IVP_H
#define _IVP_H 
#include "ANALYTIC.h"
#include <vector>
namespace iRRAM{
  template <unsigned int d, class T>
  class IVP{
      private:
        std::vector<std::shared_ptr<ANALYTIC<d,T>>> f;
        // precomputed values for the coefficients in the recursion formula
        std::vector<std::vector<T>> a;
        // ai[i,n] contains the coefficient a(i,i+n) since a(i,n)=0 for i>n
        std::vector<std::vector<std::vector<T>>> ai;
        REAL get_a(const int v,const int n, const int i);
        void init(){
          a.resize(d);
          ai.resize(d);
        };
      public:
        IVP(std::initializer_list<std::shared_ptr<ANALYTIC<d,T>>> f ): f(f) {init();};
        IVP(std::initializer_list<ANALYTIC<d,T>> fl) {
          init();
          for(auto fun : fl){
            f.push_back(std::shared_ptr<ANALYTIC<d,T>>(new ANALYTIC<d,T>(fun)));
          }
        };
        REAL get_coeff(const int v, const int n);
        REAL get_r(const int i) const {return f[i]->get_r();}
        REAL get_M(const int i) const {return f[i]->get_r()*f[i]->get_M();}
  };
  
  template<unsigned int d, class T>
  REAL IVP<d,T>::get_a(const int v, const int i, const int n){
    auto& av=ai[v];
    if(i > n) return 0;
    int m=n-i;
    if(av.size() > i && av[i].size() > m) return av[i][m];
    if(av.size() <= i){
      if(i>0)
        get_a(v,i-1, n); // make sure av.size() = i
      // add empty vector at position i
      av.push_back(std::vector<REAL>());
    }
    if(av[i].size() <= m){
      if(m > 0)
        get_a(v,i, m-1); // make sure av[i].size() = m
      av[i].push_back(0);
    }
    if(i == 0 && n==0) av[i][m] = 1;
    else if(i > 0){
      for(unsigned long j=0; j<=n;j++){
          av[i][m]+=get_a(v,i-1,n-j)*get_coeff(v,j);
      } 
    }
    return av[i][m];
  }

  template<unsigned int d, class T>
  REAL IVP<d,T>::get_coeff(const int v, const int l){
    auto& av=a[v];
    if(av.size() > l) return av[l];
    // make sure av.size() = l
    if(l > 0)
      get_coeff(v,l-1);
    av.push_back(0);
    if(l > 0)
    {
      for(int k=0; k<l;k++)
      {
        for(auto w : partitions(l-1-k, d-1)){
          for(auto I : bounded_count(w))
          {
            I.insert(I.begin(),k);
            auto c = f[v]->get_coeff(I);
            for(int i=0; i<d-1;i++)
            {
              c *= get_a(i, I[i+1], w[i]);
            }
            av[l] += c;
          }
        }
      }
      av[l] /= REAL(l);
    }
    return av[l];
  };

  template <unsigned int d, class T>
  std::vector<ANALYTIC<1,T>> solve(std::shared_ptr<IVP<d,T>> P){
    std::vector<ANALYTIC<1,T>> ans;
    for(int i=0; i<d-1; i++)
    {
      auto series = [P,i] (const std::vector<unsigned long>& v) {
        return P->get_coeff(i,(int)v[0]);
      };
      ans.push_back(ANALYTIC<1,T>(series, P->get_r(i), P->get_M(i)));
    }
    return ans;
  }
}
#endif
