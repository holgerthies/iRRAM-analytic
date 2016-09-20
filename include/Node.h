/*-------------------------------------------------
 * Abstract class for nodes in expression tree
 ------------------------------------------------*/
#ifndef NODE_H
#define NODE_H
#include <memory>
#include "POWERSERIES.h"
namespace iRRAM
{
  enum class ANALYTIC_OPERATION {ANALYTIC, ADDITION,SCALAR_ADDITION, SUBTRACTION, MULTIPLICATION,SCALAR_MULTIPLICATION, INVERSION, DERIVATIVE, COMPOSITION, IVP, SINE, COSINE, EXPONENTIATION,CONTINUATION};

  // forward declaration 
  template<class R, class... Args>
  class ANALYTIC;

  template<class R, class... Args>
  class Node
  {
  public:
    virtual ~Node() = default;
    virtual R evaluate(const Args&... args) const = 0;
    virtual std::shared_ptr<ANALYTIC<R, Args...>> to_analytic() const = 0;
    virtual ANALYTIC_OPERATION get_type() const = 0;
    template<class... IArgs>
    R evaluate_vector(const std::vector<R>& X, int pos, IArgs... iargs)
    {
      return evaluate_vector(X, pos+1, iargs..., X[pos]);
      
    }
    R evaluate_vector(const std::vector<R>& X, int pos, Args... iargs)
    {
      return evaluate(iargs...);
      
    }
    R evaluate(const std::vector<R>& X)
    {
      return evaluate_vector(X,0);
      
    }
  };

  template<class R, class... Args>
  class BinaryNode : public Node<R, Args...>
  {
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  protected:
    node_ptr lhs, rhs;
  public:
    virtual ~BinaryNode() = default;
    
    BinaryNode(const node_ptr& lhs, const node_ptr& rhs):
      lhs(lhs), rhs(rhs) 
    {
    }
  };
  
  template<class R, class... Args>
  class ScalarNode : public Node<R, Args...>
  {
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  protected:
    node_ptr node;
    R scalar;
  public:
    virtual ~ScalarNode() = default;
    
    ScalarNode(const node_ptr& node, const R& scalar):
      node(node), scalar(scalar)
    {
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
