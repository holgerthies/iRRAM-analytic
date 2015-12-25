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
        T get_a(const int v,const int n, const int i);
        void init(){
          a.resize(d);
          ai.resize(d);
        };
        T get_coeff_rec(const int v, const int l, const int k, const int upper,const int size, std::vector<unsigned long>& N);
	template <unsigned int n>
	  T get_coeff_rec2(const std::shared_ptr<POWERSERIES<n, T>>& pwr,const int pos, std::vector<unsigned long>& N);
	template <unsigned int n>
	T get_coeff_rec2(const std::shared_ptr<T>& x, const int pos, std::vector<unsigned long>& N);
    
      public:
        IVP(std::initializer_list<std::shared_ptr<ANALYTIC<d,T>>> f ): f(f) {init();};
        IVP(std::initializer_list<ANALYTIC<d,T>> fl) {
          init();
          for(auto fun : fl){
            f.push_back(std::shared_ptr<ANALYTIC<d,T>>(new ANALYTIC<d,T>(fun)));
          }
        };
        T get_coeff(const int v, const int n);
        REAL get_r(const int i) const {return f[i]->get_r();}
        REAL get_M(const int i) const {return f[i]->get_r()*f[i]->get_M();}
  };
  
  template<unsigned int d, class T>
  T IVP<d,T>::get_a(const int v, const int i, const int n){
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
    template<unsigned int n>
      T IVP<d,T>::get_coeff_rec2(const std::shared_ptr<POWERSERIES<n, T>>& pwr,const int pos, std::vector<unsigned long>& N){
      T ans=0;
      for(int i=0; i<=N[pos]; i++){
	ans += get_a(pos, i, N[pos])*get_coeff_rec2<n-1>((*pwr)[i], pos+1, N);
      }
      return ans;
    }


  template<unsigned int n, class T>
  template<unsigned int d>
    T IVP<n,T>::get_coeff_rec2(const std::shared_ptr<T>& x, const int pos, std::vector<unsigned long>& N){
    return *x;
  }


  template<unsigned int d, class T>
  T IVP<d,T>::get_coeff_rec(const int v, const int l, const int k, const int upper,const int size, std::vector<unsigned long>& N){
    using std::vector;
    if(size == d-1){
      return get_coeff_rec2<d-1>((*(f[v]->pwr))[k], 0, N);
      
    }
    T ans=0;
    int start = (size == d-2) ? upper : 0;
    for(int i=start; i<=upper; i++){
      N[size] = i;
      ans += get_coeff_rec(v, l, k, upper-i,size+1, N );
    }
    return ans;
  }

  template<unsigned int d, class T>
  T IVP<d,T>::get_coeff(const int v, const int l){
    using std::vector;
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
        vector<unsigned long> N(d-1, 0);
        av[l] += get_coeff_rec(v, l, k, l-k-1, 0, N);
      }
      av[l] /= REAL(l);
    }
    return av[l];
  }
  
  template <unsigned int d, class T>
  std::vector<ANALYTIC<1,T>> solve(std::shared_ptr<IVP<d,T>> P){
    using std::min;
    std::vector<ANALYTIC<1,T>> ans;
    for(int i=0; i<d-1; i++)
    {
      std::function<T(unsigned long)> series = [P,i] (const unsigned long n) {
        return P->get_coeff(i,n);
      };
      ans.push_back(ANALYTIC<1,T>(series, min(P->get_r(i), P->get_r(i)/P->get_M(i)), P->get_M(i)));
    }
    return ans;
  }
}
#endif
