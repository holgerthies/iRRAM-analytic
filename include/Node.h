/*-------------------------------------------------
 * Abstract class for nodes in expression tree
 ------------------------------------------------*/
#ifndef NODE_H
#define NODE_H
#include <memory>
#include "POWERSERIES.h"
#include "tutil.h"
namespace iRRAM
{
  enum class ANALYTIC_OPERATION {ANALYTIC, ADDITION,SCALAR_ADDITION, SUBTRACTION, MULTIPLICATION,SCALAR_MULTIPLICATION, INVERSION, DERIVATIVE, COMPOSITION, IVP, SINE, COSINE, EXPONENTIATION,CONTINUATION,TRANSPOSITION, POLYNOMIAL, PRUNE};

  // forward declaration 
  template<class R, class... Args>
  class ANALYTIC;
  template<class... Args>
  REAL real_max(const REAL& x, const Args&... args)
  {
    return maximum(x, real_max(args...));
  }

  template<>
  REAL real_max(const REAL& x)
  {
    return x;
  }

  template<class R, class... Args>
  class Node
  {
  private:
    using pwr_series = POWERSERIES<sizeof...(Args),R>;
    using pwr_series_ptr = std::shared_ptr<pwr_series>;
  public:
    mutable pwr_series_ptr pwr;
    Node()
    {
    }

    virtual R evaluate(const Args&... args) const {
      REAL r=real_max(args...);
      r = minimum(0.9*this->get_r(), maximum(r+1, 2*r));
      REAL M =this->get_M(r);
      return PWRSERIES_IMPL::evaluate(get_pwr(), r, M, args...);
    };
    virtual ~Node() = default;
    //virtual R evaluate(const Args&... args) const = 0;
    //virtual std::shared_ptr<ANALYTIC<R, Args...>> to_analytic() const = 0;
    virtual ANALYTIC_OPERATION get_type() const = 0;
    virtual std::string to_string() const  = 0;
    virtual REAL get_M(const REAL& r) const = 0;
    virtual REAL get_r() const = 0;
    virtual R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>&) const = 0;
    virtual size_t get_size() const {
      return 1;
    }

    template<class... IArgs>
    R evaluate_vector(const std::vector<R>& X, int pos, IArgs... iargs) const
    {
      return evaluate_vector(X, pos+1, iargs..., X[pos]);
      
    }
    R evaluate_vector(const std::vector<R>& X, int pos, Args... iargs) const
    {
      return evaluate(iargs...);
      
    }
    R evaluate(const std::vector<R>& X) const
    {
      return evaluate_vector(X,0);
      
    }
    
    virtual R get_coefficient_v(const tutil::n_tuple<sizeof...(Args),size_t>& t) const 
    {
      return get_coefficient(t);
    }

    pwr_series_ptr get_pwr() const
    {
      if(!pwr){
        auto coeff_lambda = [this] (const auto... params) {  return this->get_coefficient_v(std::make_tuple(params...));};
        auto coeff_fun= tutil::repeat_fun<sizeof...(Args), R, const size_t>(coeff_lambda);
        pwr = std::make_shared<POWERSERIES<sizeof...(Args),REAL>>(coeff_fun);
      }
      return pwr;
    }
    R get_coefficient_cached(const tutil::n_tuple<sizeof...(Args),size_t>& t) const{
      return PWRSERIES_IMPL::get_coeff<sizeof...(Args), R>(*get_pwr(), t);
    }
  };

  template<class R, class... Args>
  class BinaryNode : public Node<R, Args...>
  {
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  public:
    node_ptr lhs, rhs;
  public:
    virtual ~BinaryNode() = default;

    virtual std::string to_string() const override
    {
      return "BIN_OP("+lhs->to_string()+","+rhs->to_string()+")";
    }

    virtual size_t get_size() const override{
      return 1+lhs->get_size()+rhs->get_size();
    }
    
    BinaryNode(const node_ptr& lhs, const node_ptr& rhs):
      lhs(lhs), rhs(rhs) 
    {
    }
  };
  
  template<class R, class... Args>
  class ScalarNode : public Node<R, Args...>
  {
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  public:
    node_ptr node;
    R scalar;
  public:
    virtual ~ScalarNode() = default;
    
    ScalarNode(const node_ptr& node, const R& scalar):
      node(node), scalar(scalar)
    {
    }
    virtual std::string to_string() const override
    {
      return "SCALAR_OP("+node->to_string()+","+std::to_string(scalar.as_double())+")";
      
    }

    virtual size_t get_size() const override{
      return 1+node->get_size();
    }

  };

  template<size_t d, class T>
  class SERIES_OPERATOR
  {
  public:
    virtual ~SERIES_OPERATOR() = default;
    virtual std::shared_ptr<POWERSERIES<d-1,T>> get_coeff(const unsigned long) const = 0;
  };

  template<size_t d, class T>
  class SERIES_BINARY_OPERATOR : public SERIES_OPERATOR<d,T>
  {
    using pwrseries_ptr = std::shared_ptr<POWERSERIES<d,T>>;
  protected:
    pwrseries_ptr lhs, rhs;
  public:
    virtual ~SERIES_BINARY_OPERATOR() = default;
    
    SERIES_BINARY_OPERATOR(const pwrseries_ptr& lhs, const pwrseries_ptr& rhs):
      lhs(lhs), rhs(rhs) 
    {
    }
  };

  template<size_t d, class T>
  class SERIES_SCALAR_OPERATOR : public SERIES_OPERATOR<d,T>
  {
    using pwrseries_ptr = std::shared_ptr<POWERSERIES<d,T>>;
  protected:
    pwrseries_ptr series;
    T scalar;
  public:
    virtual ~SERIES_SCALAR_OPERATOR() = default;
    
    SERIES_SCALAR_OPERATOR(const pwrseries_ptr& series, const T& scalar):
      series(series), scalar(scalar)
    {
    }
  };
  
  template<size_t d, class T>
  std::shared_ptr<POWERSERIES<d, T>> get_series(std::shared_ptr<SERIES_OPERATOR<d,T>> op) 
  {
    return std::make_shared<POWERSERIES<d,T>>([op] (unsigned long n) {return op->get_coeff(n);});
  }
  
}
#endif
