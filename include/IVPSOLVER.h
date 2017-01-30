/*-------------------------------------------------
* Class for SERIES_IVP_SOLVER solver
 ------------------------------------------------*/
#ifndef IVP_SOLVER_H
#define IVP_SOLVER_H
#include <vector>
#include "CONTINUATION.h"
namespace iRRAM
{
  const int MAX_DAG_SIZE = 10;
  const int MAX_DAG_SIZE_LARGE = 100000;
  const int ENCLOSURE_MAX_ORDER = 35;
  // for debugging
  struct DEBUG_INFORMATION
  {
    int steps;
    int order;
  };
  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> get_next_f(const int index, const std::shared_ptr<Node<R,Args...>>& f, const std::vector<std::shared_ptr<Node<R, Args...>>>& Fs)
  {
    auto ans = pderive(f, 0,1)*Fs[0];
    for(int i=1; i<Fs.size(); i++){
      ans = ans+pderive(f, i, 1)*Fs[i];
    }
    ans = (REAL(1)/REAL(index))*ans;
    simplify(ans);
    return ans;
    
  }

  template <class R, class... Args>
  class BASIC_IVPSOLVER {
  private:
    REAL partial_sum(const int i, const REAL& x, const int m)
    {
      REAL ans = 0;
      REAL x_pow = 1;
      for(int k=0; k<m; k++){
        ans += this->get_coefficient(i, std::make_tuple<size_t>(k))*x_pow;
        x_pow *= x;
      }
      return ans;
    }

    REAL enclose(const int i, const REAL& r, const REAL& enc) const
    {
      REAL r_pow = r;
      REAL M = abs(this->get_coefficient(i, 0));
      REAL tM = M + r_pow*this->get_f(i,0)->get_M_root(enc);
      int k;
      int order = ENCLOSURE_MAX_ORDER;
      if(ACTUAL_STACK.actual_prec > -100){
        order = min(ENCLOSURE_MAX_ORDER, 5);
      }
      for(k=1; k<order; k++){
        M += abs(this->get_coefficient(i, std::make_tuple<size_t>(k)))*r_pow;
        r_pow *= r;
        REAL cM = M + r_pow*this->get_f(i,k)->get_M_root(enc);
        if(positive(cM-tM, -2)){
          break;
        }
        tM = cM;
      }
      //std::cout << "k = " << k << std::endl;
      
      //std::cout <<  this->get_f(i,order-1)->to_string()<< "\n";
      //std::cout << i << " "<<enc.as_double()<< "  " <<  M.as_double() << " " << this->get_f(i,order)->get_M_root(enc).as_double() <<std::endl;
      return tM;
    }

  protected:
    std::vector<std::shared_ptr<Node<R, Args...>>> F;
    std::vector<R> Y;
    mutable std::vector<std::vector<std::shared_ptr<Node<R, Args...>>>> Fc;
    mutable std::vector<std::vector<R>> coeffs;
    int ind;
    int enclosure_order;
    REAL q, radius;
  public:

    void set_radius()
    {
      REAL min_r = 1000000;
      // get minimum r for all functions
      int index = 0;
      REAL M=0;
      REAL new_r = 1;
      q = 0;
      for(const auto& f : F){
        if(positive(abs(Y[index])+1-f->get_r_cached(), 1)){
          std::cout << "Warning: radius of F might be too small to continue" << std::endl;
        }
        min_r = minimum(min_r, f->get_r_cached());
        q = maximum(q, abs(Y[index]));
        index++;
      }
      q = minimum(q+1, q+min_r/10);
      radius = q;
      for(int i=0; i<F.size(); i++)
      {
        REAL enc = enclose(i, radius, q);
        
        while( positive(enc-q, -1)){
          radius /= 2;
          enc = enclose(i, radius, q);
        }
      }
      //std::cout << radius.as_double() << "\n";
    }

    auto get_F(const int index) -> decltype(F[0])
    {
      return F[index];
    }

    REAL get_r() const{
      return radius;
    }
    REAL get_M(const int ind, const REAL& r) const {
      return enclose(ind, r,q );
    }
    
    virtual R get_coefficient(const int ind, const tutil::n_tuple<1,size_t>& idx) const = 0;
    virtual std::shared_ptr<Node<R, Args...>> get_f(const int ind, const int n) const = 0;

    size_t read_coefficients(const int ind){
      return coeffs[ind].size();
    }

    virtual void transpose(const std::vector<R>& Z) = 0;
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
  class IVPSOLVER : public BASIC_IVPSOLVER<R,Args...>  {
  public:
    IVPSOLVER(const std::vector<std::shared_ptr<Node<R,Args...>>>& F, const std::vector<R>& Y)
    {
      this->F = F;
      this->Y = Y;
      this->coeffs = std::vector<std::vector<R>>(F.size(), std::vector<R>());
      for(const auto& f : F){
        this->Fc.push_back(std::vector<std::shared_ptr<Node<R,Args...>>>{f});
      }
      this->set_radius();
    }

    std::shared_ptr<Node<R, Args...>> get_f(const int ind, const int n) const override
    {
      while(this->Fc[ind].size() <= n){
        int i = this->Fc[ind].size();
        //std::cout << "("<<ind << ", " << i<<") " << std::endl;
        this->Fc[ind].push_back(get_next_f(i+1, this->Fc[ind][i-1], this->F));
        if(this->Fc[ind][i]->get_size() > MAX_DAG_SIZE_LARGE)
          this->Fc[ind][i] = prune(this->Fc[ind][i]);
        //std::cout << this->Fc[ind][n]->get_size() <<" "<<this->Fc[ind][n]->node_number()<< "("<<ind << ", " << i<<") " << std::endl;
      }
      return this->Fc[ind][n];
    }

    R get_coefficient(const int ind, const tutil::n_tuple<1,size_t>& idx) const override
    {
      int n = std::get<0>(idx);
      if(this->coeffs[ind].size() == 0)
        this->coeffs[ind].push_back(this->Y[ind]);

      while(this->coeffs[ind].size() <= n){
        REAL m=this->get_f(ind, this->coeffs[ind].size()-1)->evaluate(this->Y);
        this->coeffs[ind].push_back(m);
        
      }
      return this->coeffs[ind][n];
    }

    void transpose(const std::vector<R>& Z) override
    {
      this->Y = Z;
      for(auto& v : this->coeffs){
        v.clear();
      }
      
      this->set_radius();
    }
  };

  template <class R, class... Args>
  class IVPSOLVER_REC : public BASIC_IVPSOLVER<R,Args...>  {
  private:
    std::vector<std::shared_ptr<Node<R,Args...>>> F_orig;
  public:
    IVPSOLVER_REC(const std::vector<std::shared_ptr<Node<R,Args...>>>& F, const std::vector<R>& Y)
    {
      this->F_orig = F;
      this->coeffs = std::vector<std::vector<R>>(F.size(), std::vector<R>());
      this->transpose(Y);
    }

    std::shared_ptr<Node<R, Args...>> get_f(const int ind, const int n) const override
    {
      while(this->Fc[ind].size() <= n){
        int i = this->Fc[ind].size();
        this->Fc[ind].push_back(get_next_f(i+1, this->Fc[ind][i-1], this->F));
        if(this->Fc[ind][i]->get_size() > MAX_DAG_SIZE)
          this->Fc[ind][i] = prune(this->Fc[ind][i]);
        
        //std::cout << this->Fc[ind][i]->node_number() << "("<<ind << ", " << i<<") " << std::endl;
      }
      return this->Fc[ind][n];
    }

    R get_coefficient(const int ind, const tutil::n_tuple<1,size_t>& idx) const override
    {
      int n = std::get<0>(idx);
      if(this->coeffs[ind].size() == 0)
        this->coeffs[ind].push_back(this->Y[ind]);

      while(this->coeffs[ind].size() <= n){
        REAL m=this->get_f(ind, this->coeffs[ind].size()-1)->get_zero_coefficient();
        this->coeffs[ind].push_back(m);
        
      }
      return this->coeffs[ind][n];
    }

    void transpose(const std::vector<R>& Z) override
    {
      this->Y = std::vector<R>(Z.size(), 0);
      for(auto& v : this->coeffs){
        v.clear();
      }
      this->Fc.clear();
      this->F.clear();
      for(const auto& f : this->F_orig){
        auto ft = iRRAM::transpose(f, Z);
        simplify(ft);
        this->F.push_back(ft);
        this->Fc.push_back(std::vector<std::shared_ptr<Node<R,Args...>>>{ft});
      }
      this->set_radius();
    }
  };
  template <class R, class... Args>
  class IVP : public Node<R, R>{
  private:
    std::shared_ptr<BASIC_IVPSOLVER<R,Args...>> solver;
    int ind;
  public:
    IVP(const std::shared_ptr<BASIC_IVPSOLVER<R,Args...>>& solver, const int ind):
      solver(solver), ind(ind)
    {
    }

    REAL get_r() const override {
      return solver->get_r();
    }
    REAL get_M(const REAL& r) const override {
      return solver->get_M(ind, r);
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
    R t0;
    std::vector<R> y;
  };
  

  // solve ODE system by taylor method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> ivp_solve_cs(const IVPSYSTEM<R, Args...>& S, const R& max_time, const bool output, const int solver_type, DEBUG_INFORMATION& debug)
  {
    single_valued code;
    std::vector<R> Y(S.y);
    std::vector<R> Y0(S.y);
    std::vector<R> Z(S.y.size(), 0);
    int iter=0, order=0;
    std::shared_ptr<BASIC_IVPSOLVER<R, Args...>> solver;
    if(solver_type == 0){
      // make sure Y(0)=0 in the beginning
      decltype(S.F) Fc;
      for(const auto& f: S.F){
        auto fc = iRRAM::transpose(f, Y0);
        simplify(fc);
        Fc.push_back(fc);
      }
      solver = std::make_shared<IVPSOLVER<R,Args...>>(Fc, Y);
    }
    else{
      Z = Y0;
      solver = std::make_shared<IVPSOLVER_REC<R,Args...>>(S.F, Y);
    }
    REAL t = S.t0;
    while(positive(max_time-t, -3)){
      iter++;
      solver->transpose(Z);
      REAL h;
      for(int i=0; i<S.F.size(); i++){
        auto F = std::make_shared<IVP<R, Args...>>(solver, i);
        int d = 5;//-ACTUAL_STACK.actual_prec*S.F.size()/10;
        h=F->get_r_cached()/d;
        //REAL p = minimum(h, max_time-Y[0]);
        if(solver_type == 0)
        {
          Z[i] = F->evaluate_root(h);
          Y[i] = Y0[i]+Z[i];
        }
        else
        {
          Z[i] += F->evaluate_root(h);
          Y[i] = Z[i];
        }
        order = max(order, std::dynamic_pointer_cast<IVP<R,Args...>>(F)->read_coefficients());
      }
      t += h;
      if(output){
        iRRAM::cout <<iter << " " << t;
        for(auto y : Y)
          iRRAM::cout << " " << y;
        iRRAM::cout << std::endl;
        
      }
    }
    debug.steps = iter;
    debug.order = order;
    return  Y;
  }

  // solve ODE system by taylor method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> ivp_solve_co(const IVPSYSTEM<R, Args...>& S, const R& max_time,const int order, const bool output, const int solver_type, DEBUG_INFORMATION& debug)
  {
    single_valued code;
    std::vector<R> Y(S.y);
    std::vector<R> Y0(S.y);
    std::vector<R> Z(S.y.size(), 0);
    int iter=0;
    std::shared_ptr<BASIC_IVPSOLVER<R, Args...>> solver;
    if(solver_type == 0){
      // make sure Y(0)=0 in the beginning
      decltype(S.F) Fc;
      for(const auto& f: S.F){
        auto fc = iRRAM::transpose(f, Y0);
        simplify(fc);
        Fc.push_back(fc);
      }
      solver = std::make_shared<IVPSOLVER<R,Args...>>(Fc, Y);
    }
    else{
      solver = std::make_shared<IVPSOLVER_REC<R,Args...>>(S.F, Y);
    }
    REAL t = S.t0;
    while(positive(max_time-t, -3)){
      iter++;
      solver->transpose(Y);
      REAL h=1;
      for(int i=0; i<S.F.size(); i++){
          auto F = std::make_shared<IVP<R, Args...>>(solver, i);
          h = minimum(h, F->get_r());
      }
      auto trunc_error = real_to_error(1);
      while (trunc_error.exponent >= ACTUAL_STACK.actual_prec){
        h /= 2;
        sizetype_exact(trunc_error);
        for(int i=0; i<S.F.size(); i++){
          auto F = std::make_shared<IVP<R, Args...>>(solver, i);
          REAL error = abs(solver->get_f(i,order)->get_M_root(abs(Z[i])+1))*power(h, order);
          if(solver_type == 0)
          {
            Z[i] = F->evaluate_root(h);
          }
          else
          {
            Z[i] = Y[i]+F->evaluate_root(h);
          }
          auto local_trunc_error = real_to_error(error);
          sizetype sum_error, local_error;
          Z[i].geterror(sum_error);
          sizetype_add(local_error, sum_error, local_trunc_error);
          Z[i].seterror(local_error);
          
          sizetype_max(trunc_error,trunc_error, local_trunc_error);
        }
      }
      for(int i=0; i<S.F.size(); i++){
        Y[i] = Z[i];
      }
      
      t += h;
      if(output){
        iRRAM::cout <<iter << " " << t;
        for(auto y : Y)
          iRRAM::cout << " " << y;
        iRRAM::cout << std::endl;
        
      }
    }
    debug.steps = iter;
    debug.order = order;
    return  Y;
  }


  // solve ODE system by taylor method and evaluate at some point x
  template<class R, class... Args>
  std::vector<R> ivp_solve_mixed(const IVPSYSTEM<R, Args...>& S, const R& max_time, const bool output, const int solver_type, DEBUG_INFORMATION& debug)
  {
    single_valued code;
    std::vector<R> Y(S.y);
    std::vector<R> Y0(S.y);
    std::vector<R> Z(S.y.size(), 0);
    int iter=0, order=0;
    std::shared_ptr<BASIC_IVPSOLVER<R, Args...>> solver;
    if(solver_type == 0){
      // make sure Y(0)=0 in the beginning
      decltype(S.F) Fc;
      for(const auto& f: S.F){
        auto fc = iRRAM::transpose(f, Y0);
        simplify(fc);
        Fc.push_back(fc);
      }
      solver = std::make_shared<IVPSOLVER<R,Args...>>(Fc, Y);
    }
    else{
      Z = Y0;
      solver = std::make_shared<IVPSOLVER_REC<R,Args...>>(S.F, Y);
    }
    REAL t = S.t0;
    while(positive(max_time-t, -3)){
      iter++;
      solver->transpose(Z);
      REAL h;
      for(int i=0; i<S.F.size(); i++){
        auto F = std::make_shared<IVP<R, Args...>>(solver, i);
        int d = -ACTUAL_STACK.actual_prec*S.F.size()/16;
        h=F->get_r_cached()/d;
        //REAL p = minimum(h, max_time-Y[0]);
        if(solver_type == 0)
        {
          Z[i] = F->evaluate_root(h);
          Y[i] = Y0[i]+Z[i];
        }
        else
        {
          Z[i] += F->evaluate_root(h);
          Y[i] = Z[i];
        }
        order = max(order, std::dynamic_pointer_cast<IVP<R,Args...>>(F)->read_coefficients());
      }
      t += h;
      if(output){
        iRRAM::cout <<iter << " " << t;
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
