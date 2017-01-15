/*-------------------------------------------------
* Class for SERIES_IVP_SOLVER solver
 ------------------------------------------------*/
#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H
#include <vector>
#include "CONTINUATION.h"
namespace iRRAM
{
  const int MAX_DAG_SIZE = 1;
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
  class BASIC_IVPSOLVER {
  protected:
    std::vector<std::shared_ptr<Node<R, Args...>>> F;
    std::vector<R> Y;
    mutable std::vector<std::vector<std::shared_ptr<Node<R, Args...>>>> Fc;
    mutable std::vector<std::vector<R>> coeffs;
    int ind;
    REAL first_radius, radius;
  public:

    void set_radius()
    {
      REAL min_r = 1;
      // get minimum r for all functions
      for(const auto& f : F){
        min_r = minimum(min_r, f->get_r_cached());
      }

      REAL new_r = 0;
      radius = -1;
      while(positive(new_r-radius,-5)){
        REAL M =0;
        int i=0;
        for(const auto& f : F){
          i++;
          M = maximum(M, f->get_M(min_r+abs(Y[0])));
        }
        radius = new_r;
        new_r = min_r/M;
        //std::cout << new_r.as_double() << "> " << min_r.as_double() << " < " << M.as_double() << std::endl;
        min_r /= 2;
      }
    }

    auto get_F(const int index) -> decltype(F[0])
    {
      return F[index];
    }

    REAL get_r() const{
      return radius;
    }
    REAL get_M(const int ind, const REAL& r) const {
      return r*F[ind]->get_M(r+abs(Y[0]))+abs(Y[ind+1]);
    }
    
    virtual R get_coefficient(const int ind, const tutil::n_tuple<1,size_t>& idx) const = 0;

    size_t read_coefficients(const int ind){
      return coeffs[ind].size();
    }

    virtual size_t get_size() const {
      size_t ans=1;
      for(auto& f : F)
        ans+=f->get_size();
      return ans;
    }
    void reset_visited() const 
    {
      for(auto& f : F){
        f->reset_visited();
      }
    }

    int count_nodes() const 
    {
      int n=0;
      for(auto& f : F){
        n += f->count_nodes();
      }
      return n;
    }

    std::string to_string(const int ind) const 
    {
      std::string ans= "IVP(";
      for(const auto& f : F){
        ans += f->to_string()+", ";
      }
      ans += "i = "+std::to_string(ind)+")";
      return ans;
    }
  };

  template <class R, class... Args>
  class IVPSOLVER_REC : public BASIC_IVPSOLVER<R,Args...>  {
  private:
    std::vector<std::shared_ptr<Node<R, Args...>>> F_original;
  public:
    IVPSOLVER_REC(const std::vector<std::shared_ptr<Node<R,Args...>>>& F, const std::vector<R>& Y)
    {
      this->F_original = F;
      // Y0 is zero
      this->Y = std::vector<R>(Y.size(), R());
      this->coeffs = std::vector<std::vector<R>>(F.size(), std::vector<R>());
      transpose(Y);
    }
    R get_coefficient(const int ind, const tutil::n_tuple<1,size_t>& idx) const override
    {
      int n = std::get<0>(idx);
      if(this->coeffs[ind].size() == 0)
        this->coeffs[ind].push_back(this->Y[ind+1]);

      while(this->Fc[ind].size() <= n){
        int i = this->Fc[ind].size();
        this->Fc[ind].push_back(get_next_f(i+1, this->Fc[ind][i-1], this->F));
        if(this->Fc[ind][i]->get_size() > MAX_DAG_SIZE)
          this->Fc[ind][i] = prune(this->Fc[ind][i]);
      }
      while(this->coeffs[ind].size() <= n){
        
        REAL m=this->Fc[ind][this->coeffs[ind].size()-1]->evaluate(this->Y);
        this->coeffs[ind].push_back(m);
        
      }
      return this->coeffs[ind][n];
    }

    void transpose(const std::vector<R>& Z)
    {
      this->Y = Z;
      for(auto& v : this->coeffs){
        v.clear();
      }
      this->Fc.clear();
      this->F.clear();
      for(const auto& f : F_original){
        auto ft = iRRAM::transpose(f, Z);
        this->F.push_back(ft);
        this->Fc.push_back(std::vector<std::shared_ptr<Node<R,Args...>>>{ft});
      }
      this->set_radius();
    }
  };

  template <class R, class... Args>
  class IVPSOLVER : public BASIC_IVPSOLVER<R,Args...>  {
  public:
    IVPSOLVER(const std::vector<std::shared_ptr<Node<R,Args...>>>& F, const std::vector<R>& Y)
    {
      this->F = F;
      this->Y = Y;
      this->coeffs = std::vector<std::vector<R>>(F.size(), std::vector<R>());
      this->set_radius();
      for(const auto& f : F){
        this->Fc.push_back(std::vector<std::shared_ptr<Node<R,Args...>>>{f});
      }
    }
    R get_coefficient(const int ind, const tutil::n_tuple<1,size_t>& idx) const override
    {
      int n = std::get<0>(idx);
      if(this->coeffs[ind].size() == 0)
        this->coeffs[ind].push_back(this->Y[ind+1]);

      while(this->Fc[ind].size() <= n){
        int i = this->Fc[ind].size();
        this->Fc[ind].push_back(get_next_f(i+1, this->Fc[ind][i-1], this->F));
        if(this->Fc[ind][i]->get_size() > MAX_DAG_SIZE)
          this->Fc[ind][i] = prune(this->Fc[ind][i]);
        std::cout << "("<<ind << ", " << i<<") " << std::endl;
      }
      while(this->coeffs[ind].size() <= n){
        
        REAL m=this->Fc[ind][this->coeffs[ind].size()-1]->evaluate(this->Y);
        this->coeffs[ind].push_back(m);
        
      }
      return this->coeffs[ind][n];
    }

    void transpose(const std::vector<R>& Z)
    {
      this->Y = Z;
      for(auto& v : this->coeffs){
        v.clear();
      }
      
      this->set_radius();
    }
  };

  template <class R, class... Args>
  class IVP : public Node<R, R>{
  private:
    REAL r_improved, r_orig;
    std::shared_ptr<BASIC_IVPSOLVER<R,Args...>> solver;
    int ind;
  public:
    IVP(const std::shared_ptr<BASIC_IVPSOLVER<R,Args...>>& solver, const int ind):
      solver(solver), ind(ind)
    {
      r_improved = solver->get_r();
      r_orig = r_improved;
    }

    REAL get_r() const override {
      return r_improved;
    }
    REAL get_M(const REAL& r) const override {
      REAL rp = minimum(r, r_orig);
      return solver->get_M(ind, rp);
    }
    
    void improve_radius()
    {
      REAL q = 0.55*r_improved;
      REAL M = get_M(r_orig);
      // upper bound for F(q)
      REAL ub = (abs(this->approximate(20,q)))+1;
      //std::cout << M.as_double() << " " << ub.as_double() << std::endl;
      REAL rp = q;
      auto F = solver->get_F(ind);
      while(positive(M- (ub+rp*F->get_M(rp+q)) ,-5)){
        r_improved = q+rp;
        rp *= 2;
        
      }
    }

    virtual size_t get_size() const override{
      return solver->get_size();
    }
    void reset_visited() const override
    {
      if(this->visited){
        solver->reset_visited();
        this->visited = false;
      }
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        return 1+solver->count_nodes();
      }
      return 0;
    }

    R get_coefficient(const tutil::n_tuple<1,size_t>& idx) const override
    {
      return solver->get_coefficient(ind, idx);
    }

    size_t read_coefficients(){
      return solver->read_coefficients(ind);
    }

    std::string to_string() const override
    {
      return solver->to_string(ind);
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
  std::shared_ptr<Node<R,R>> ivp_solve(const std::vector<std::shared_ptr<Node<R,Args...>>>& F,const std::vector<R>& Y, const int i)
  {
    auto solver = std::make_shared<IVPSOLVER<R,Args...>>(F, Y);
    return std::make_shared<IVP<R, Args...>>(solver, i);
  }

  // solve ODE system by taylor method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> solve_taylor(const IVPSYSTEM<R, Args...>& S, const R& max_time, bool output,  DEBUG_INFORMATION& debug)
  {
    if(!output)
      single_valued code;
    std::vector<R> Y(S.y);
    int iter=0, order=0;

    auto solver = std::make_shared<IVPSOLVER<R,Args...>>(S.F, Y);
    while(Y[0] < max_time){
      iter++;
      solver->transpose(Y);
      REAL h;
      for(int i=0; i<S.F.size(); i++){
        auto F = std::make_shared<IVP<R, Args...>>(solver, i);
        // std::cout << "r: " << F->get_r().as_double() <<std::endl;
        // F->improve_radius();
        // std::cout << "r improved: " << F->get_r().as_double() <<std::endl;
        // F->improve_radius();
        // std::cout << "r improved: " << F->get_r().as_double() <<std::endl;
        int d =-ACTUAL_STACK.actual_prec*S.F.size()/10;
        h=F->get_r_cached()/d;
        //REAL p = minimum(h, max_time-Y[0]);
        Y[i+1] = F->evaluate_root(h);
        order = max(order, std::dynamic_pointer_cast<IVP<R,Args...>>(F)->read_coefficients());
      }
      Y[0] += h;
      if(output){
        iRRAM::cout <<iter;
        for(auto y : Y)
          iRRAM::cout << " " << y;
        iRRAM::cout << std::endl;
        
      }
    }
    debug.steps = iter;
    debug.order = order;
    return  Y;
  }



} // namespace iRRAM


#endif
