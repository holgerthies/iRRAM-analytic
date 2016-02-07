/*-------------------------------------------------
* Class for division
 ------------------------------------------------*/
#ifndef DIVISION_H
#define DIVISION_H
#include "ANALYTIC.h"
namespace iRRAM
{
  template <size_t d, class T>
  class SERIES_INVERSION : public SERIES_OPERATOR<d,T>
  {
  private:
    mutable std::vector<std::shared_ptr<POWERSERIES<d-1,T>>> cache;
    std::shared_ptr<POWERSERIES<d,T>> series;
  public:
    SERIES_INVERSION(const std::shared_ptr<POWERSERIES<d,T>>& series):
      series(series) 
    {
    }
    
    std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long n) const override 
    {
      auto sz = cache.size();
      if(n >= sz){
        cache.resize(n+1);
      }
      for(auto j=sz; j<=n; j++){
        if(j == 0){
          std::shared_ptr<SERIES_OPERATOR<d-1,T>> inverter = std::make_shared<SERIES_INVERSION<d-1, T>>((*this->series)[0]);
          cache[0] = get_series(inverter);
        }
        else{
	  auto sp = std::make_shared<POWERSERIES<d-1,T>>(T());
	  cache[j] = sp;
	  for(int k=0; k<=j; k++){
            std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplier = std::make_shared<SERIES_MULTIPLICATION<d-1,T>>((*this->series)[k],cache[j-k] );
            std::shared_ptr<SERIES_OPERATOR<d-1,T>> adder = std::make_shared<SERIES_ADDITION<d-1,T>>(cache[j],get_series(multiplier) );
	    cache[j] = get_series(adder);
	  }
          std::shared_ptr<SERIES_OPERATOR<d-1,T>> multiplier = std::make_shared<SERIES_MULTIPLICATION<d-1,T>>(cache[0], cache[j]);
	  std::shared_ptr<SERIES_OPERATOR<d-1,T>> scalar_multiplier = std::make_shared<SERIES_SCALAR_MULTIPLICATION<d-1,T>>(get_series(multiplier), T(-1));
	  cache[j] = get_series(scalar_multiplier);
        }
      }
      return cache[n];
    }

  };

  template<class T>
  class SERIES_INVERSION<1,T> : public SERIES_OPERATOR<1,T>
  {
  private:
    mutable std::vector<T> cache;
    std::shared_ptr<POWERSERIES<1,T>> series;
  public:
    SERIES_INVERSION(const std::shared_ptr<POWERSERIES<1,T>>& series):
      series(series) 
    {
    }
    
    std::shared_ptr<T> get_coeff(const unsigned long n) const override 
    {
      auto sz = cache.size();
      if(n >= sz){
        cache.resize(n+1);
      }
      for(auto j=sz; j<=n; j++){
        if(j == 0){
          cache[0] = T(1)/series->get(0);
        }
        else{
	  cache[j] = T();
	  for(int k=0; k<=j; k++){
            cache[j] += series->get(k)*cache[j-k];
	  }
	  cache[j] = -cache[0]*cache[j];
        }
      }
      return std::make_shared<T>(cache[n]);
    }
  };
  

  template <class R, class... Args>
  class INVERSION : public Node<R, Args...>{
  private:
    std::shared_ptr<Node<R,Args...>> node;
  public:
    INVERSION(const std::shared_ptr<Node<R, Args...>>& node):
      node(node) 
    {
    }
    R evaluate(const Args&... args) const override;
    std::shared_ptr<ANALYTIC<R,Args...>> to_analytic() const override;
    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::INVERSION;
    }
  };

  // member definitions 
  template <class R, class... Args>
  R INVERSION<R, Args...>::evaluate(const Args&... args) const
  {
    return 1/this->node->evaluate(args...);
  }

  template <class R, class... Args>
  std::shared_ptr<ANALYTIC<R,Args...>> INVERSION<R,Args...>::to_analytic() const
    {
      auto f = this->node->to_analytic();
      std::shared_ptr<SERIES_OPERATOR<sizeof...(Args), R>> division= std::make_shared<SERIES_INVERSION<sizeof...(Args), R>>(f->get_series());
      auto inv_pwr = get_series(division);
      auto newr = f->get_r();
      // find r' such that the function does not have a zero for all z \in B_r'
      REAL lowerbound=-1;// lower bound for the function
      auto f0 = constant_coefficient(f->get_series());
      while(!positive(lowerbound, -10)){
        newr /= 2;
        lowerbound = abs(f0)-int(sizeof...(Args))*f->get_M()*newr/f->get_r();
      }
      auto new_M = R(1)/lowerbound;
      return std::make_shared<ANALYTIC<R, Args...>> (inv_pwr, new_M, newr);
    }

  // division operators
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

  

} // namespace iRRAM


#endif
