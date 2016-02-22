/*-------------------------------------------------
* Class for SERIES_IVP_SOLVER solver
 ------------------------------------------------*/
#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H
#include "ANALYTIC.h"
namespace iRRAM
{

  template <size_t d, class T>
  class SERIES_IVP_SOLVER
  {
  private:
    std::vector<std::shared_ptr<POWERSERIES<d,T>>> series;
    // precomputed values for the coefficients in the recursion formula
    std::vector<std::vector<T>> a;
    // ai[i,n] contains the coefficient a(i,i+n) since a(i,n)=0 for i>n
    std::vector<std::vector<std::vector<T>>> ai;
    T get_a(const int v,const int n, const int i);
    T get_coeff_rec(const int v, const int l, const int k, const int upper,const int size, std::vector<unsigned long>& N);
    template <size_t n>
    T get_coeff_rec2(const std::shared_ptr<POWERSERIES<n, T>>& pwr,const int pos, std::vector<unsigned long>& N);
    template <size_t n>
    T get_coeff_rec2(const std::shared_ptr<T>& x, const int pos, std::vector<unsigned long>& N);
  public:
    SERIES_IVP_SOLVER(const std::vector<std::shared_ptr<POWERSERIES<d,T>>>& series ): series(series) {
      a.resize(d);
      ai.resize(d);
      
    };
    T get_coeff(const int v, const int n);
  };

  template<size_t d, class T>
  T SERIES_IVP_SOLVER<d,T>::get_a(const int v, const int i, const int n){
    auto& av=ai[v];
    if(i > n) return 0;
    int m=n-i;
    if(av.size() > i && av[i].size() > m) return av[i][m];
    if(av.size() <= i){
      if(i>0)
        get_a(v,i-1, n); // make sure av.size() = i
      // add empty vector at position i
      av.push_back(std::vector<T>());
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

  template<size_t d, class T>
  template<size_t n>
  T SERIES_IVP_SOLVER<d,T>::get_coeff_rec2(const std::shared_ptr<POWERSERIES<n, T>>& pwr,const int pos, std::vector<unsigned long>& N){
    T ans=0;
    for(int i=0; i<=N[pos]; i++){
      ans += get_a(pos, i, N[pos])*get_coeff_rec2<n-1>((*pwr)[i], pos+1, N);
    }
    return ans;
  }


  template<size_t n, class T>
  template<size_t d>
  T SERIES_IVP_SOLVER<n,T>::get_coeff_rec2(const std::shared_ptr<T>& x, const int pos, std::vector<unsigned long>& N){
    return *x;
  }


  template<size_t d, class T>
  T SERIES_IVP_SOLVER<d,T>::get_coeff_rec(const int v, const int l, const int k, const int upper,const int size, std::vector<unsigned long>& N){
    using std::vector;
    if(size == d-1){
      return get_coeff_rec2<d-1>((*(series[v]))[k], 0, N);
      
    }
    T ans=0;
    int start = (size == d-2) ? upper : 0;
    for(int i=start; i<=upper; i++){
      N[size] = i;
      ans += get_coeff_rec(v, l, k, upper-i,size+1, N );
    }
    return ans;
  }

  template<size_t d, class T>
  T SERIES_IVP_SOLVER<d,T>::get_coeff(const int v, const int l){
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

  template<size_t d, class T>
  std::shared_ptr<POWERSERIES<1, T>> get_series(const int v, std::shared_ptr<SERIES_IVP_SOLVER<d,T>> op) 
  {
    auto coeff_fun = std::function<T(const unsigned long)>([v,op] (const unsigned long n)  {return op->get_coeff(v, n);});
    return std::make_shared<POWERSERIES<1,T>>(coeff_fun);
  }

  template <class R, class... Args>
  class IVP_SOLVER
  {
  private:
    std::shared_ptr<SERIES_IVP_SOLVER<sizeof...(Args),R>> series;
    std::vector<std::shared_ptr<ANALYTIC<R, R>>> solutions;
    REAL M,r;
  public:
    IVP_SOLVER(std::initializer_list<std::shared_ptr<Node<R, Args...>>> funs)
    {
      
      M = (*funs.begin())->to_analytic()->get_M();
      r = (*funs.begin())->to_analytic()->get_r();
      std::vector<std::shared_ptr<POWERSERIES<sizeof...(Args),R>>> pwrs;
      for(auto p : funs){
        
        auto f = p->to_analytic();
        M = maximum(M, f->get_M());
        r = minimum(r, f->get_r());
        pwrs.push_back(f->get_series());
      }
      auto sol_M = r;
      auto sol_r = minimum(r, r/M);
      
      series = std::make_shared<SERIES_IVP_SOLVER<sizeof...(Args), R>>(pwrs);
      for(int i=0; i<sizeof...(Args)-1; i++){
        solutions.push_back(std::make_shared<ANALYTIC<R, R>>(get_series(i, series), sol_M, sol_r));
      }
    }
    R evaluate(const int v, const R& arg) const
    {
      return solutions[v]->evaluate(arg);
    }
    std::shared_ptr<ANALYTIC<R, R>> to_analytic(const int v)
    {
      return solutions[v];
    }
  };
  
  template <class R, class... Args>
  class IVP : public Node<R, R>{
    using solver_ptr = std::shared_ptr<IVP_SOLVER<R, Args...>>;
  private:
    solver_ptr solver;
    int v;
  public:
    IVP(const solver_ptr& solver, const int v):
      solver(solver), v(v)
    {
    }

    R evaluate(const R& arg) const override{
      return solver->evaluate(v, arg);
    };
    std::shared_ptr<ANALYTIC<R,R>> to_analytic() const override{
      return solver->to_analytic(v);
    };

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::IVP;
    }
  };

  template <class R, class... Args>
  std::vector<std::shared_ptr<Node<R,R>>> ivp_solve(std::initializer_list<std::shared_ptr<Node<R,Args...>>> F)
  {
    auto solver = std::make_shared<IVP_SOLVER<R,Args...>>(F);
    std::vector<std::shared_ptr<Node<R,R>>> ans;
    for(int i=0; i<sizeof...(Args)-1; i++){
      
      std::shared_ptr<Node<R,R>> sol = std::make_shared<IVP<R,Args...>>(solver, i);
      ans.push_back(sol);
    }
    return ans;
  }

} // namespace iRRAM


#endif
