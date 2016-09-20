/*-------------------------------------------------
* Class for SERIES_IVP_SOLVER solver
 ------------------------------------------------*/
#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H
#include "ANALYTIC.h"
#include <vector>
namespace iRRAM
{
  // for debugging
  struct DEBUG_INFORMATION
  {
    int steps;
    int order;
  };
  
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
    IVP_SOLVER(std::vector<std::shared_ptr<Node<R, Args...>>> funs)
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

  
  template <class... Args>
  REAL abs_maximum(const REAL& x, const Args&... xs )
  {
    return maximum(abs(x), abs_maximum(xs...));
  }

  REAL abs_maximum(const REAL& x)
  {
    return abs(x);
  }

  // class to specify IVP system together with initial values
  template <class R, class... Args>
  class IVPSYSTEM
  {
  public:
    std::vector<std::shared_ptr<Node<R, Args...>>> F;
    std::vector<R> y;
  };
  

  template <class R, class... Args>
  std::vector<std::shared_ptr<Node<R,R>>> ivp_solve(std::vector<std::shared_ptr<Node<R,Args...>>> F)
  {
    auto solver = std::make_shared<IVP_SOLVER<R,Args...>>(F);
    std::vector<std::shared_ptr<Node<R,R>>> ans;
    for(int i=0; i<sizeof...(Args)-1; i++){
      std::shared_ptr<Node<R,R>> sol = std::make_shared<IVP<R,Args...>>(solver, i);
      ans.push_back(sol);
    }
    return ans;
  }

  // solve ODE system by taylor method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> solve_taylor(const IVPSYSTEM<R, Args...>& S, const R& max_time,  DEBUG_INFORMATION& debug)
  {
    std::vector<R> Y(S.y);
    // transposed functions
    std::vector<std::shared_ptr<Node<R, Args...>>> Fc(S.F);
    REAL M = S.F[0]->to_analytic()->get_M();
    REAL r0 = S.F[0]->to_analytic()->get_r();
    REAL r = r0;
    int iter=0;
    while(Y[0] < max_time){
      iter++;
      
      for(int i=0; i<Fc.size(); i++){
        // transpose function such that y(0)=0
        Fc[i] = continuation(S.F[i], M,r, Y );
      }
      
      auto F = ivp_solve(Fc);
      REAL h=F[0]->to_analytic()->get_r()/8;
      //iRRAM::cout << Y[0] << std::endl;
      REAL max_y = Y[0];
      for(int i=0; i<Fc.size(); i++){
        REAL p = minimum(h, max_time-Y[0]);
        Y[i+1] += F[i]->evaluate(p);
        max_y = maximum(max_y, Y[i+1]);
      }

      r = r0-max_y;
      Y[0] += h;
      
      //iRRAM::cout << Y[0]<<"  "<<Y[2]<<std::endl;
    }
    debug.steps = iter;
    debug.order = S.F[0]->to_analytic()->get_series()->known_coeffs();
    
    return  Y;
  }

  template<class R, class... Args>
  bool heun_step(const std::vector<std::shared_ptr<Node<R, Args...>>>& F, const REAL& max_time, std::vector<R>& Y)
  {
    REAL M = F[0]->to_analytic()->get_M();
    REAL r = F[0]->to_analytic()->get_r();
    for(const auto& f : F){
      M = maximum(M, f->to_analytic()->get_M());
      r = minimum(r, f->to_analytic()->get_r());
    }
  
    sizetype eval_error, trunc_error;
    sizetype_exact(eval_error);
    REAL error=1;
    REAL error_factor = (10*M+3*M*M)/(3*r*r);
    REAL h = r/8-Y[0];
    bool stop=false;
    if(!positive(max_time-Y[0]-h, -3)){
      h = max_time-Y[0];
      stop=true;
    }
    
    
    std::vector<REAL> yp(Y.size());
    std::vector<REAL> ny(Y.size());
    bool halfed=false;
    
    do{
      if(halfed) stop=false;
      yp[0] = Y[0]+h;
      ny[0] = Y[0]+h;
      for(int i=0; i<F.size(); i++){
        yp[i+1] = Y[i+1]+h*F[i]->evaluate(Y);
      }
      for(int i=0; i<F.size(); i++){
        ny[i+1] = (Y[i+1]+yp[i+1]+h*F[i]->evaluate(yp))/2;
        sizetype eval_error0, eval_error1;
        sizetype_copy(eval_error1, eval_error);
        ny[i+1].geterror(eval_error0);
        sizetype_max(eval_error, eval_error0, eval_error1);
      }
      error = power(h,3)*error_factor;
      sizetype error_error, error_vsize;
      error.geterror(error_error);
      error.getsize(error_vsize);
      sizetype_add(trunc_error,error_vsize,error_error);
      h /= 2;
      halfed = true;
      
      
    } while ( 
      (trunc_error.exponent >= ACTUAL_STACK.actual_prec ) );
    
    
    Y[0] = ny[0];
    
    for(int i=1; i<Y.size(); i++){
      sizetype local_error, total_error;
      Y[i] = ny[i];
      Y[i].geterror(local_error);
      sizetype_add(total_error, local_error, trunc_error);
      Y[i].seterror(total_error);
      
    }
    return stop;
  }

  
  template<class R, class... Args>
  bool heun_step_aff(const std::vector<std::shared_ptr<Node<R, Args...>>>& F, const REAL& max_time, std::vector<AAREAL>& Y)
  {
    REAL M = F[0]->to_analytic()->get_M();
    REAL r = F[0]->to_analytic()->get_r();
    for(const auto& f : F){
      M = maximum(M, f->to_analytic()->get_M());
      r = minimum(r, f->to_analytic()->get_r());
    }
  
    sizetype eval_error, trunc_error;
    sizetype_exact(eval_error);
    REAL error=1;
    REAL error_factor = (10*M+3*M*M)/(3*r*r);
    AAREAL h = r/8-Y[0];
    bool stop=false;
    if(positive(max_time-Y[0].to_real()-h.to_real(), -3)){
      h = max_time-Y[0];
      stop=true;
    }
    
    
    std::vector<AAREAL> yp(Y.size());
    std::vector<REAL> ypr(Y.size());
    std::vector<AAREAL> ny(Y.size());
    std::vector<REAL> Yr(Y.size());

    for(int i = 0; i < Y.size(); i++){
      Yr[i] = Y[i].to_real();
    }

    bool halfed=false;
    
    do{
      if(halfed) stop=false;
      yp[0] = Y[0]+h;
      ypr[0] = yp[0].to_real();
      ny[0] = Y[0]+h;
      for(int i=0; i<F.size(); i++){
        yp[i+1] = Y[i+1]+h*F[i]->evaluate(Yr);
        ypr[i+1] = yp[i+1].to_real();
      }
      for(int i=0; i<F.size(); i++){
        ny[i+1] = (Y[i+1]+yp[i+1]+h*F[i]->evaluate(ypr))*(REAL(1)/2);
        sizetype eval_error0, eval_error1;
        sizetype_copy(eval_error1, eval_error);
        ny[i+1].to_real().geterror(eval_error0);
        sizetype_max(eval_error, eval_error0, eval_error1);
      }
      error = power(h.to_real(),3)*error_factor;
      sizetype error_error, error_vsize;
      error.geterror(error_error);
      error.getsize(error_vsize);
      sizetype_add(trunc_error,error_vsize,error_error);
      h = REAL(1)/REAL(2)*h;
      
      halfed = true;
      
      
    } while ( 
      (trunc_error.exponent >= ACTUAL_STACK.actual_prec ) );
    
    
    Y[0] = ny[0];
    
    for(int i=1; i<Y.size(); i++){
      Y[i] = ny[i];
      Y[i].add_error(error);
      
    }
    return stop;
  }

  template<class R, class... Args>
  bool euler_step(const std::vector<std::shared_ptr<Node<R, Args...>>>& F, const REAL& max_time, std::vector<R>& Y)
  {
    REAL M = F[0]->to_analytic()->get_M();
    REAL r = F[0]->to_analytic()->get_r();
    for(const auto& f : F){
      M = maximum(M, f->to_analytic()->get_M());
      r = minimum(r, f->to_analytic()->get_r());
    }
  
    sizetype eval_error, trunc_error;
    sizetype_exact(eval_error);
    REAL error_factor=M*M;
    REAL error=1;
    REAL h = r/8-Y[0]+REAL(0.001);
    bool stop=false;
    if(!positive(max_time-Y[0]-h, -3)){
      h = max_time-Y[0];
      stop=true;
    }
    
    
    std::vector<REAL> ny(Y.size());
    
    bool halfed=false;
    
    do{
      if(halfed) stop=false;
      ny[0] = Y[0]+h;
      for(int i=0; i<F.size(); i++){
        ny[i+1] = Y[i+1]+h*F[i]->evaluate(Y);
        sizetype eval_error0, eval_error1;
        sizetype_copy(eval_error1, eval_error);
        ny[i+1].geterror(eval_error0);
        sizetype_max(eval_error, eval_error0, eval_error1);
      }
      error = power(h,2)*error_factor/(r-Y[0]-2*h);
      trunc_error = real_to_error(error);
      h /= 2;
      halfed=true;
    } while ( 
      (trunc_error.exponent >= ACTUAL_STACK.actual_prec ) );
    
    
    Y[0] = ny[0];
    
    for(int i=1; i<Y.size(); i++){
      sizetype local_error, total_error;
      Y[i] = ny[i];
      Y[i].geterror(local_error);
      sizetype_add(total_error, local_error, trunc_error);
      Y[i].seterror(total_error);
    }
    return stop;
  }

  
  // solve ODE system by heun's method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> solve_euler(const IVPSYSTEM<R, Args...>& S, const R& max_time,  DEBUG_INFORMATION& debug)
  {
    std::vector<R> Y(S.y);
    
    int iter=0;
    REAL before = Y[0];
    
    while(!euler_step(S.F, max_time, Y)){
      //std::cout << iter << " " << Y[0].as_double()  << "\n";
      ++iter;
    }
    debug.steps = iter;
    debug.order = 1;
    return Y;
  }

  
  // solve ODE system by heun's method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> solve_heun(const IVPSYSTEM<R, Args...>& S, const R& max_time,  DEBUG_INFORMATION& debug)
  {
    std::vector<REAL> Y(S.y);
    int iter=0;
    
    while(!heun_step(S.F, max_time, Y)){
      ++iter;
      //std::cout <<iter << " " << Y[0].as_double() << "\n";
      // std::cout << "error:" << Y[1].error.mantissa << "&"<<Y[1].error.exponent << std::endl;
      
    }
    debug.steps = iter;
    debug.order = 2;
    
    return Y;
  }

  template<class R, class... Args>
  std::vector<R> solve_heun_aff(const IVPSYSTEM<R, Args...>& S, const R& max_time,  DEBUG_INFORMATION& debug)
  {
    std::vector<AAREAL> Y(S.y.size());
    std::vector<REAL> Yr(S.y);
    for(int i=0; i<Y.size(); i++){
      Y[i] = S.y[i];
      
    }
    int iter=0;
    
    while(!heun_step_aff(S.F, max_time, Y)){
      ++iter;
      if(iter % 5 == 0)
        for(auto& y : Y) y.clean();
      
      //std::cout <<iter << " " << Y[0].to_real().as_double() << "\n";
      // std::cout << "error:" << Y[1].error.mantissa << "&"<<Y[1].error.exponent << std::endl;
      
    }
    debug.steps = iter;
    debug.order = 2;
    for(int i=0; i<Y.size();i++){
      Yr[i] = Y[i].to_real();
    }
    std::cout << "returning" << std::endl;
    
    return Yr;
  }

  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> get_next_f(const int index, const std::shared_ptr<Node<R,Args...>>& f, const std::vector<std::shared_ptr<Node<R, Args...>>>& Fs)
  {
    auto ans = pderive(f, 0,1);
    for(int i=0; i<Fs.size(); i++){
            ans = ans+pderive(f, i+1, 1)*Fs[i];
    }
    ans = (REAL(1)/REAL(index))*ans;
    return ans;
    
  }

  template<class R, class... Args>
  bool taylor_step(const int order, const std::vector<std::shared_ptr<Node<R, Args...>>>& F, const REAL& max_time, std::vector<R>& Y)
  {
    std::vector<std::shared_ptr<Node<R, Args...>>> Fc(F);
    std::vector<std::shared_ptr<Node<R, Args...>>> Fco(F);

    std::vector<REAL> Z(Y.size());
    REAL M = F[0]->to_analytic()->get_M();
    REAL r = F[0]->to_analytic()->get_r();
    for(const auto& f : F){
      M = maximum(M, f->to_analytic()->get_M());
      r = minimum(r, f->to_analytic()->get_r());
    }
  
    for(int i=0; i<Fc.size(); i++){
        // transpose function such that y(0)=0
        Fc[i] = continuation(F[i], M,r, Y );
        Fco[i] = continuation(F[i], M,r, Y );
    }
    sizetype eval_error, trunc_error;
    sizetype_exact(eval_error);
    REAL error_factor=power(2,order+1)*M/power(r,order+1);
    REAL error=1;
    REAL h = 0.1;
    do{
      h /=2;
      error = power(h,order+1)*error_factor;
      
      trunc_error = real_to_error(error);
    } while ( 
      (trunc_error.exponent >= ACTUAL_STACK.actual_prec ) );
  
    
    bool stop=false;
    
    std::vector<REAL> ny(Y);
    
    
    
    ny[0] = Y[0]+h;
    REAL hpow = 1;
    h = minimum(max_time-Y[0],h);
    
    for(int i=1; i<=order; i++){
      hpow *= h;
      for(int j=0; j<F.size(); j++){
        REAL m=Fc[j]->evaluate(Z);
        Fc[j] = get_next_f(i+1, Fc[j], Fco);
        REAL t=Fc[j]->evaluate(Z);
        
        ny[j+1] += m*hpow;
      }
    }

    Y[0] = ny[0];
    
    for(int i=1; i<Y.size(); i++){
      sizetype local_error, total_error;
      Y[i] = ny[i];
      Y[i].geterror(local_error);
      sizetype_add(total_error, local_error, trunc_error);
      Y[i].seterror(total_error);
    }
    if(Y[0] > max_time) return true;
    
    return stop;
  }

  template<class R, class... Args>
  std::vector<R> solve_taylor_deriv(const IVPSYSTEM<R, Args...>& S, const R& max_time, DEBUG_INFORMATION& debug)
  {
    std::vector<REAL> Y(S.y);
    int iter=0;
    int order=20;
    while(!taylor_step(order, S.F, max_time, Y)){
      ++iter;
      //std::cout <<iter << " " << Y[0].as_double() <<" " << Y[1].as_double()<< "\n";
      
      // std::cout << "error:" << Y[1].error.mantissa << "&"<<Y[1].error.exponent << std::endl;
      
    }
    debug.steps = iter;
    debug.order = order;
    
    return Y;
  }



} // namespace iRRAM


#endif
