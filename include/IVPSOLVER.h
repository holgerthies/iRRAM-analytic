/*-------------------------------------------------
* Class for SERIES_IVP_SOLVER solver
 ------------------------------------------------*/
#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H
#include <vector>
namespace iRRAM
{
  const int MAX_DAG_SIZE = 10;
  // for debugging
  struct DEBUG_INFORMATION
  {
    int steps;
    int order;
  };
  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> get_next_f(const int index, const std::shared_ptr<Node<R,Args...>>& f, const std::vector<std::shared_ptr<Node<R, Args...>>>& Fs)
  {
    auto ans = pderive(f, 0,1);
    for(int i=0; i<Fs.size(); i++){
      ans = ans+pderive(f, i+1, 1)*Fs[i];
    }
    ans = (REAL(1)/REAL(index))*ans;
    simplify(ans);
    return ans;
    
  }

  template <class R, class... Args>
  class IVP : public Node<R, R>{
  private:
    std::vector<std::shared_ptr<Node<R, Args...>>> F;
    mutable std::shared_ptr<Node<R, Args...>> Fc;
    mutable std::vector<R> coeffs;
    int ind;
    REAL radius;
  public:
    IVP(const std::vector<std::shared_ptr<Node<R,Args...>>>& F, const int ind):
      F(F),Fc(F[ind]), ind(ind)
    {
      radius = 1;
      REAL M =0;
      for(const auto& f : F){
        radius = minimum(radius, f->get_r());
      }
      for(const auto& f : F){
        M = maximum(M, f->get_M(radius));
      }
      radius = radius/M;
    }

    REAL get_r() const override {
      return radius;
    }
    REAL get_M(const REAL& r) const override {
      return r*F[ind]->get_M(r);
    }
    
    virtual size_t get_size() const override{
      size_t ans=1;
      for(auto& f : F)
        ans+=f->get_size();
      return ans;
    }

    R get_coefficient(const tutil::n_tuple<1,size_t>& idx) const override
    {
      int n = std::get<0>(idx);
      tutil::n_tuple<sizeof...(Args),size_t> Zt;
      if(coeffs.size() == 0)
        coeffs.push_back(0);
      while(coeffs.size() <= n){
        REAL m=Fc->get_coefficient(Zt);
        Fc = get_next_f(coeffs.size()+1, Fc, F);
        if(Fc->get_size() > MAX_DAG_SIZE)
          Fc = prune(Fc);
        coeffs.push_back(m);
      }
      return coeffs[n];
    }

    size_t read_coefficients(){
      return coeffs.size();
    }

    std::string to_string() const override
    {
      std::string ans= "IVP(";
      for(const auto& f : F){
        ans += f->to_string()+", ";
      }
      ans += "i = "+std::to_string(ind)+")";
      return ans;
    }
    
    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::IVP;
    }
  };

  
  // class to specify IVP system together with initial values
  template <class R, class... Args>
  class IVPSYSTEM
  {
  public:
    std::vector<std::shared_ptr<Node<R, Args...>>> F;
    std::vector<R> y;
  };
  

  template <class R, class... Args>
  std::shared_ptr<Node<R,R>> ivp_solve(const std::vector<std::shared_ptr<Node<R,Args...>>>& F, const int i)
  {
    return std::make_shared<IVP<R, Args...>>(F, i);
  }

  // solve ODE system by taylor method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> solve_taylor(const IVPSYSTEM<R, Args...>& S, const R& max_time,  DEBUG_INFORMATION& debug)
  {
    single_valued code;
    std::vector<R> Y(S.y);
    // transposed functions
    int iter=0, order=0;
    std::vector<std::shared_ptr<Node<R, Args...>>> Fc(S.F.size());
    while(Y[0] < max_time){
      iter++;
      for(int i=0; i<Fc.size(); i++){
        // transpose function such that y(0)=0
        Fc[i] = transpose(S.F[i], Y );
        simplify(Fc[i]);
      }

      REAL h;
      for(int i=0; i<Fc.size(); i++){
        auto F = ivp_solve(Fc,i);
        int d =-ACTUAL_STACK.actual_prec/2;
        h=F->get_r()/d;
        REAL max_y = Y[0];
        REAL p = minimum(h, max_time-Y[0]);
        Y[i+1] += F->evaluate(p);
        order = max(order, std::dynamic_pointer_cast<IVP<R,Args...>>(F)->read_coefficients());
      }
      Y[0] += h;

      //std::cout <<iter << " " << Y[0].as_double() <<" " << Y[1].as_double()<< "\n";
    }
    debug.steps = iter;
    debug.order = order;
    return  Y;
  }


  template<class R, class... Args>
  bool taylor_step(const int order, const std::vector<std::shared_ptr<Node<R, Args...>>>& F, const REAL& max_time, std::vector<R>& Y)
  {
    std::vector<std::shared_ptr<Node<R, Args...>>> Fc(F);
    std::vector<std::shared_ptr<Node<R, Args...>>> Fco(F);

    // zero vector
    std::vector<R> Z(sizeof...(Args));
    tutil::n_tuple<sizeof...(Args),size_t> Zt;
    REAL radius = 1;
    REAL M = F[0]->get_M(1);
    for(const auto& f : F){
      M = maximum(M, f->get_M(1));
    }
  

    for(int i=0; i<Fc.size(); i++){
      // transpose function such that y(0)=0
      Fc[i] = transpose(F[i], Y );
      simplify(Fc[i]);
      Fco[i] = transpose(F[i], Y );
      simplify(Fco[i]);
    }
    REAL error;
    sizetype trunc_error;
    bool stop=false;
    
    std::vector<REAL> ny(Y);
    
    auto coefficients = std::vector<std::vector<R>>(F.size(), std::vector<R>(order+1));
    for(int i=1; i<=order; i++){
      //std::cout << i << std::endl;
      for(int j=0; j<F.size(); j++){
        REAL m=Fc[j]->get_coefficient(Zt);
        // std::cout << "getting coeff end" << std::endl;
        //if(i % 5 == 0) Fc[j] = prune(Fc[j]);
         // std::cout << Fc[j]->to_string() << std::endl;
         //  std::cout << std::endl;
         //  std::cout << std::endl;
        Fc[j] = get_next_f(i+1, Fc[j], Fco);
        // std::cout << std::endl;
        coefficients[j][i] = m;
      }
    }
    REAL h = minimum(0.99*radius/M, 0.99*max_time);
    do{
      error = 0;
      h /=2;
      for(int j=0; j<F.size(); j++){
        //std::cout << Fc[j]->get_r().as_double() << "::" << Fc[j]->get_M(h).as_double()<< "\n";
        error = maximum(error, Fc[j]->get_M(h)*power(h,order+1));
      }
      trunc_error = real_to_error(error);
    } while ( 
      (trunc_error.exponent >= ACTUAL_STACK.actual_prec ) );
    ny[0] = Y[0]+h;
    h = minimum(max_time-Y[0],h);

    REAL hpow = 1;

    for(int i=1; i<=order; i++){
      error = 0;
      hpow *= h;
      for(int j=0; j<F.size(); j++){
        ny[j+1] += coefficients[j][i]*hpow;
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
    
    //if(Y[0] > max_time) return true;
    
    return stop;
  }

  template<class R, class... Args>
  std::vector<R> solve_taylor_deriv(const IVPSYSTEM<R, Args...>& S, const R& max_time, DEBUG_INFORMATION& debug)
  {
    single_valued code;
    std::vector<REAL> Y(S.y);
    int iter=0;
    int order=max(60,-ACTUAL_STACK.actual_prec/3);
    while(Y[0] < max_time){
      taylor_step(order, S.F, max_time, Y);
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
