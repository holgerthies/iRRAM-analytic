/*-------------------------------------------------
* Class for division
 ------------------------------------------------*/
#ifndef DIVISION_H
#define DIVISION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  // template <size_t d, class T>
  // class SERIES_INVERSION : public SERIES_OPERATOR<d,T>
  // {
  // private:
  //   mutable std::vector<std::shared_ptr<POWERSERIES<d-1,T>>> cache;
  //   std::shared_ptr<POWERSERIES<d,T>> series;
  // public:
  //   SERIES_INVERSION(const std::shared_ptr<POWERSERIES<d,T>>& series):
  //     series(series) 
  //   {
  //   }
    
  //   std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
  //   {
  //     auto sz = cache.size();
  //     if(n >= sz){
  //       cache.resize(n+1);
  //     }
  //     for(auto j=sz; j<=n; j++){
  //       if(j == 0){
  //         std::shared_ptr<SERIES_OPERATOR<d-1,T>> inverter = std::make_shared<SERIES_INVERSION<d-1, T>>((*this->series)[0]);
  //         cache[0] = get_series(inverter);
  //       }
  //       else{
  //         auto sp = std::make_shared<POWERSERIES<d-1,T>>(T());
  //         cache[j] = sp;
  //         for(int k=0; k<=j; k++){
  //           std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplier = std::make_shared<SERIES_MULTIPLICATION<d-1,T>>((*this->series)[k],cache[j-k] );
  //           std::shared_ptr<SERIES_OPERATOR<d-1,T>> adder = std::make_shared<SERIES_ADDITION<d-1,T>>(cache[j],get_series(multiplier) );
  //           cache[j] = get_series(adder);
  //         }
  //         std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplier = std::make_shared<SERIES_MULTIPLICATION<d-1,T>>(cache[0], cache[j]);
  //         std::shared_ptr<SERIES_OPERATOR<d-1,T>> scalar_multiplier = std::make_shared<SERIES_SCALAR_MULTIPLICATION<d-1,T>>(get_series(multiplier), T(-1));
  //         cache[j] = get_series(scalar_multiplier);
  //       }
  //     }
  //     return cache[n];
  //   }

  // };

  // template<class T>
  // class SERIES_INVERSION<1,T> : public SERIES_OPERATOR<1,T>
  // {
  // private:
  //   mutable std::vector<T> cache;
  //   std::shared_ptr<POWERSERIES<1,T>> series;
  // public:
  //   SERIES_INVERSION(const std::shared_ptr<POWERSERIES<1,T>>& series):
  //     series(series) 
  //   {
  //   }
    
  //   std::shared_ptr<T> get_coeff(const unsigned long n) const override 
  //   {
  //     auto sz = cache.size();
  //     if(n >= sz){
  //       cache.resize(n+1);
  //     }
  //     for(auto j=sz; j<=n; j++){
  //       if(j == 0){
  //         cache[0] = T(1)/series->get(0);
  //       }
  //       else{
  //         cache[j] = T();
  //         for(int k=0; k<=j; k++){
  //           cache[j] += series->get(k)*cache[j-k];
  //         }
  //         cache[j] = -cache[0]*cache[j];
  //       }
  //     }
  //     return std::make_shared<T>(cache[n]);
  //   }
  // };
  
  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> invert(const std::shared_ptr<Node<R,Args...>>& node);
  
  template <class R, class... Args>
  class INVERSION : public Node<R, Args...>{
  private:
    mutable std::shared_ptr<Node<R,Args...>> grad;
    mutable REAL r,f0;
  public:
    std::shared_ptr<Node<R,Args...>> node;
    INVERSION(const std::shared_ptr<Node<R, Args...>>& node):
      node(node) 
    {
    }

    R evaluate(const Args&... args) const override{
      return 1/node->evaluate_cached(args...);
    }

    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      size_t n = std::get<0>(idx);
      if(n == 0){
        auto a0 = fix_first(this->node, 0);
        return invert(a0)->get_coefficient_cached(tutil::tail(idx));
      }
      R fact=this->get_coefficient_cached(tutil::tuple_replace<0>(idx, 0));
      R bn = node->get_coefficient_cached(tutil::tuple_replace<0>(idx, 1))*this->get_coefficient_cached(tutil::tuple_replace<0>(idx, n-1));
      for(int i=2; i<=n; i++){
        bn = bn+node->get_coefficient_cached(tutil::tuple_replace<0>(idx, i))*this->get_coefficient_cached(tutil::tuple_replace<0>(idx, n-i));
      }
      return REAL(-1)*fact*bn;
    }

    REAL get_r() const override {
      if(!grad){
       // find r' such that the function does not have a zero for all
       // z \in B_r'
        grad = pderive(node, 0, 1);
        for(int i=1; i<sizeof...(Args);i++)
          grad = grad+pderive(node, i, 1);
        simplify(grad);
        r=node->get_r_cached();
        std::vector<R> Z(sizeof...(Args));
        f0 = node->evaluate(Z);
        REAL lowerbound=-1;// lower bound for the function
        while(!positive(lowerbound, -10)){
          r /= 2;
          lowerbound = abs(f0)-grad->get_M_root(r)*r;
        }
      } 
      return r;
    }

    void reset_visited() const override
    {
      if(this->visited){
        this->visited = false;
        node->reset_visited();
      }
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        int n=1+node->count_nodes();
        return n;
      }
      return 0;
    }

    std::string to_string() const override
    {
      return "1/("+this->node->to_string()+")";
      
    }

    REAL get_M(const REAL& r) const override {
      get_r();
      REAL lowerbound = abs(f0)-grad->get_M_cached(r)*r;
      return 1/lowerbound;
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::INVERSION;
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
    }
  };

  template <class R>
  class INVERSION<R,R> : public Node<R, R>{
  private:
    mutable std::shared_ptr<Node<R,R>> df;
    mutable REAL r,f0;
  public:
    std::shared_ptr<Node<R,R>> node;
    INVERSION(const std::shared_ptr<Node<R, R>>& node):
      node(node) 
    {
    }

    R evaluate(const R& x) const override{
      return 1/node->evaluate_cached(x);
    }

    R get_coefficient(const std::tuple<size_t>& idx) const override
    {
      
      size_t n = std::get<0>(idx);
      if(n == 0){
        auto a0 = this->node->get_coefficient_cached(idx);
        return REAL(1)/a0;
      }
      R bn=0;
      for(int i=1; i<=n; i++){
        bn = bn+node->get_coefficient_cached(std::make_tuple(i))*this->get_coefficient_cached(std::make_tuple(n-i));
      }
      return REAL(-1)*this->get_coefficient_cached(std::make_tuple(0))*bn;
    }

    REAL get_r() const override {
      if(!df){
       // find r' such that the function does not have a zero for all
       // z \in B_r'
        df = pderive(node, 0, 1);
        simplify(df);
        r=node->get_r_cached();
        f0 = node->evaluate(0);
        REAL lowerbound=-1;// lower bound for the function
        while(!positive(lowerbound, -10)){
          r /= 2;
          lowerbound = abs(f0)-df->get_M_root(r)*r;
        }
      } 
      return r;
    }

    void reset_visited() const override
    {
      if(this->visited){
        this->visited = false;
        node->reset_visited();
      }
    }

    int count_nodes() const override
    {
      if(!this->visited){
        this->visited = true;
        int n=1+node->count_nodes();
        return n;
      }
      return 0;
    }

    std::string to_string() const override
    {
      return "1/("+this->node->to_string()+")";
      
    }

    REAL get_M(const REAL& r) const override {
      get_r();
      REAL lowerbound = abs(f0)-df->get_M_cached(r)*r;
      return 1/lowerbound;
    }

    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::INVERSION;
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
    }
  };

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator/(const std::shared_ptr<Node<R,Args...>>& lhs,const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    
    std::shared_ptr<Node<R,Args...>> inv = std::make_shared<INVERSION<R, Args...>>(rhs);
    return inv*lhs;
  }


  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator/(const std::shared_ptr<Node<R,Args...>>& lhs,const R& rhs)
  {
    
    return std::make_shared<SCALAR_MULTIPLICATION<R, Args...>>(lhs, R(1)/rhs);
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> operator/(const R& lhs, const std::shared_ptr<Node<R,Args...>>& rhs)
  {
    std::shared_ptr<Node<R,Args...>> inv = std::make_shared<INVERSION<R, Args...>>(rhs);
    return inv*lhs;
  }

  template <class R, class... Args>
  std::shared_ptr<Node<R, Args...>> invert(const std::shared_ptr<Node<R,Args...>>& node)
  {
    std::shared_ptr<Node<R,Args...>> inv = std::make_shared<INVERSION<R, Args...>>(node);
    return inv;
  }
  

} // namespace iRRAM


#endif
