#ifndef SIMPLIFICATION_H
#define SIMPLIFICATION_H
namespace iRRAM{


  // Prune node says the everything below should not be further
  // simplified
  template <class R, class... Args>
  class PRUNE : public Node<R, Args...>{
    using node_ptr = std::shared_ptr<Node<R, Args...>>;
  public:
    node_ptr node;
    PRUNE(const node_ptr& node):
     node(node)
    {
    }

    REAL get_r() const override
    {
      return node->get_r_cached();
    }

    REAL get_M(const REAL& r) const override
    {
      return node->get_M_cached(r);
    }
    
    R get_coefficient(const tutil::n_tuple<sizeof...(Args),size_t>& idx) const override
    {
      return node->get_coefficient_cached(idx);
    }

    std::string to_string() const override
    {
      std::string ans="PRUNED_NODE("+node->to_string()+")";
      return ans;
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
    ANALYTIC_OPERATION get_type() const override
    {
      return ANALYTIC_OPERATION::PRUNE;
    }
  };

  template<class R, class... Args>
  bool simplify_step(std::shared_ptr<Node<R, Args...>>&);
  
  template<size_t v, class R, class... Args>
  struct simplify_product_rule_helper
  {
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
    {
      auto orders_t = f->orders_t;
      const int order = std::get<v>(orders_t);
      if(order == 0)
        return simplify_product_rule_helper<v+1,R,Args...>::get(f);
      auto child = std::dynamic_pointer_cast<MULTIPLICATION<R,Args...>>(f->node);
      auto ans = pderive(child->lhs, v, order)*child->rhs;
      for(int i=1; i<=order; i++){
        ans = ans+REAL(choose(order,i))*pderive(child->lhs, v, order-i)*pderive(child->rhs, v, i);
      }
      return tderive(ans, tutil::tuple_replace<v>(orders_t, 0));
    }
  };

  template<class R, class... Args>
  struct simplify_product_rule_helper<sizeof...(Args),R,Args...>
  {
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
    {
      return f->node;
    }
      
  };
    
  
  template<size_t v, class R, class... Args>
  struct simplify_chain_rule_helper
  {
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
    {
      auto orders_t = f->orders_t;
      const int order = std::get<v>(orders_t);
      if(order == 0)
        return simplify_chain_rule_helper<v+1,R,Args...>::get(f);
      auto child = std::dynamic_pointer_cast<COMPOSITION<R,Args...>>(f->node);
      auto ans = compose(pderive(child->lhs, 0, 1), child->rhs)*pderive(child->rhs, v, 1);
      return tderive(ans, tutil::tuple_replace<v>(orders_t, order-1));
    }
  };

  template<class R, class... Args>
  struct simplify_chain_rule_helper<sizeof...(Args),R,Args...>
  {
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
    {
      return f->node;
    }
  };

  template<size_t v, class R, class... Args>
  struct simplify_derivative_inverse_helper
  {
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
    {
      auto orders_t = f->orders_t;
      const int order = std::get<v>(orders_t);
      if(order == 0)
        return simplify_derivative_inverse_helper<v+1,R,Args...>::get(f);
      auto child = std::dynamic_pointer_cast<INVERSION<R,Args...>>(f->node);
      auto ans = (REAL(-1)/(child->node * child->node))*pderive(child->node, v, 1);
      return tderive(ans, tutil::tuple_replace<v>(orders_t, order-1));
    }
  };

  template<class R, class... Args>
  struct simplify_derivative_inverse_helper<sizeof...(Args),R,Args...>
  {
    static std::shared_ptr<Node<R,Args...>> get(const std::shared_ptr<DERIVATIVE<R,Args...>>& f)
    {
      return f->node;
    }
  };
  template<class R, class... Args>
  bool simplify_addition(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<ADDITION<R,Args...>>(f);
    if(node->lhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL && node->rhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto lchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->lhs);
      auto rchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->rhs);
      f = std::make_shared<poly_impl::POLY<R, Args...>>(*lchild + *rchild);
      return true;
    }
    if(node->lhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto lchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->lhs);
      //is zero?
      if(get_max_degree(*lchild) == 0)
      {
        f = node->rhs;
        return true;
      }
      // is constant?
      if(get_max_degree(*lchild) == 1)
      {
        tutil::n_tuple<sizeof...(Args), size_t> Z;
        f = node->rhs+lchild->get_coefficient_cached(Z);
        return true;
      }
    }
    if(node->rhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto rchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->rhs);
      //is zero?
      if(get_max_degree(*rchild) == 0)
      {
        f = node->lhs;
        return true;
      }
      // is constant?
      if(get_max_degree(*rchild) == 1)
      {
        tutil::n_tuple<sizeof...(Args), size_t> Z;
        f = node->lhs+rchild->get_coefficient_cached(Z);
        return true;
      }
    }
    bool lhs_changed = simplify_step(node->lhs);
    bool rhs_changed = simplify_step(node->rhs);
    return lhs_changed || rhs_changed;
  }

  template<class R, class... Args>
  bool simplify_subtraction(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<SUBTRACTION<R,Args...>>(f);
    if(node->lhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL && node->rhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto lchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->lhs);
      auto rchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->rhs);
      f = std::make_shared<poly_impl::POLY<R, Args...>>(*lchild + REAL(-1)* (*rchild));
      return true;
    }
    bool lhs_changed = simplify_step(node->lhs);
    bool rhs_changed = simplify_step(node->rhs);
    return lhs_changed || rhs_changed;
  }

  template<class R, class... Args>
  bool simplify_scalar_addition(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<SCALAR_ADDITION<R,Args...>>(f);
    switch(node->node->get_type()){
    case ANALYTIC_OPERATION::SCALAR_ADDITION:
      {
        auto child = std::dynamic_pointer_cast<SCALAR_ADDITION<R,Args...>>(node->node);
        f = std::make_shared<SCALAR_ADDITION<R,Args...>>(child->node, node->scalar+child->scalar);
        return true;
      }
    case ANALYTIC_OPERATION::POLYNOMIAL:
      {
        auto child = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->node);
        auto new_node = std::make_shared<poly_impl::POLY<R,Args...>>(node->scalar+(*child));
        f = new_node;
        return true;
      }
    default:
      return simplify_step(node->node);
    }
  }

  template<class R, class... Args>
  bool simplify_multiplication(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<MULTIPLICATION<R,Args...>>(f);
    if(node->lhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL && node->rhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto lchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->lhs);
      auto rchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->rhs);
      f = std::make_shared<poly_impl::POLY<R, Args...>>(*lchild * *rchild);
      return true;
    }
    if(node->lhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto lchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->lhs);
      //is zero?
      if(get_max_degree(*lchild) == 0)
      {
        f = std::make_shared<poly_impl::POLY<R, Args...>>();
        return true;
      }
      // is constant?
      if(get_max_degree(*lchild) == 1)
      {
        tutil::n_tuple<sizeof...(Args), size_t> Z;
        f = node->rhs*lchild->get_coefficient_cached(Z);
        return true;
      }
    }
    if(node->rhs->get_type() == ANALYTIC_OPERATION::POLYNOMIAL){
      auto rchild = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->rhs);
      //is zero?
      if(get_max_degree(*rchild) == 0)
      {
        f = std::make_shared<poly_impl::POLY<R, Args...>>();
        return true;
      }
      // is constant?
      if(get_max_degree(*rchild) == 1)
      {
        tutil::n_tuple<sizeof...(Args), size_t> Z;
        f = node->lhs*rchild->get_coefficient_cached(Z);
        return true;
      }
    }
    bool lhs_changed = simplify_step(node->lhs);
    bool rhs_changed = simplify_step(node->rhs);
    return lhs_changed || rhs_changed;
  };

  template<class R, class... Args>
  bool simplify_scalar_multiplication(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<SCALAR_MULTIPLICATION<R,Args...>>(f);
    switch(node->node->get_type()){
    case ANALYTIC_OPERATION::SCALAR_MULTIPLICATION:
      {

        auto child = std::dynamic_pointer_cast<SCALAR_MULTIPLICATION<R,Args...>>(node->node);
        auto new_node = std::make_shared<SCALAR_MULTIPLICATION<R,Args...>>(child->node, node->scalar*child->scalar);
        f = new_node;
        return true;
      }
    case ANALYTIC_OPERATION::POLYNOMIAL:
      {
        auto child = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->node);
        auto new_node = std::make_shared<poly_impl::POLY<R,Args...>>(node->scalar*(*child));
        f = new_node;
        return true;
      }
    default:
      return simplify_step(node->node);
    }
  }

  template<class R, class... Args>
  bool simplify_composition(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<COMPOSITION<R,Args...>>(f);
    return simplify_step(node->lhs) || simplify_step(node->rhs);
  }

  template<class R, class... Args>
  bool simplify_derivative_gen(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<DERIVATIVE<R,Args...>>(f);
    auto orders_t = node->orders_t;
    switch(node->node->get_type()){
    case ANALYTIC_OPERATION::ADDITION:
      {
        auto child = std::dynamic_pointer_cast<ADDITION<R,Args...>>(node->node);
        auto new_node = std::make_shared<ADDITION<R,Args...>>(tderive(child->lhs,orders_t), tderive(child->rhs, orders_t));
        f = new_node;
        return true;
      }
    case ANALYTIC_OPERATION::SUBTRACTION:
      {
        auto child = std::dynamic_pointer_cast<SUBTRACTION<R,Args...>>(node->node);
        auto new_node = std::make_shared<SUBTRACTION<R,Args...>>(tderive(child->lhs,orders_t), tderive(child->rhs, orders_t));
        f = new_node;
        return true;
      }
    case ANALYTIC_OPERATION::SCALAR_ADDITION:
      {
        auto child = std::dynamic_pointer_cast<SCALAR_ADDITION<R,Args...>>(node->node);
        auto new_node = tderive(child->node,orders_t);
        f = new_node;
        return true;
      }
    case ANALYTIC_OPERATION::MULTIPLICATION:
      {
        f = simplify_product_rule_helper<0, R, Args...>::get(std::make_shared<DERIVATIVE<R,Args...>>(*node));
        return true;
      }

    case ANALYTIC_OPERATION::SCALAR_MULTIPLICATION:
      {
        auto child = std::dynamic_pointer_cast<SCALAR_MULTIPLICATION<R,Args...>>(node->node);
        auto new_node = std::make_shared<SCALAR_MULTIPLICATION<R,Args...>>(tderive(child->node,orders_t), child->scalar);
        f = new_node;
        return true;
      } 

    case ANALYTIC_OPERATION::COMPOSITION:
      {
        f = simplify_chain_rule_helper<0, R, Args...>::get(std::make_shared<DERIVATIVE<R,Args...>>(*node));
        return true;
      }
    case ANALYTIC_OPERATION::INVERSION:
      {
        f = simplify_derivative_inverse_helper<0, R, Args...>::get(std::make_shared<DERIVATIVE<R,Args...>>(*node));
        return true;
      }
    case ANALYTIC_OPERATION::DERIVATIVE:
    {
      auto child = std::dynamic_pointer_cast<DERIVATIVE<R, Args...>>(node->node);
      auto new_node = tderive(child->node, tutil::sum_tuples(child->orders_t, node->orders_t));
      f = new_node;
      return true;
    }
    case ANALYTIC_OPERATION::POLYNOMIAL:{
      auto child = std::dynamic_pointer_cast<poly_impl::POLY<R, Args...>>(node->node);
      auto new_node = std::make_shared<poly_impl::POLY<R, Args...>>(derive(*child, node->orders_t));
      f = new_node;
      return true;
    }
    default:
      return simplify_step(node->node);
    }

  }

  template<class R, class... Args>
  bool simplify_derivative(std::shared_ptr<Node<R, Args...>>& f)
  {
    return simplify_derivative_gen(f);
  }

  template<>
  bool simplify_derivative<REAL, REAL>(std::shared_ptr<Node<REAL, REAL>>& f)
  {
    auto node = std::dynamic_pointer_cast<DERIVATIVE<REAL, REAL>>(f);
    switch(node->node->get_type()){
    case ANALYTIC_OPERATION::SINE:
      {
      auto child = std::dynamic_pointer_cast<SINE>(node->node);
      auto order = std::get<0>(node->orders_t);
      if(order % 4 == 0)
        f = sine_function;
      if(order % 4 == 1)
        f = cosine_function;
      if(order % 4 == 2)
        f = REAL(-1)*sine_function;
      if(order % 4 == 3)
        f = REAL(-1)*cosine_function;
      return true;
      }
    case ANALYTIC_OPERATION::COSINE:
      {
      auto child = std::dynamic_pointer_cast<SINE>(node->node);
      auto order = std::get<0>(node->orders_t);
      if(order % 4 == 0)
        f = cosine_function;
      if(order % 4 == 1)
        f = REAL(-1)*sine_function;
      if(order % 4 == 2)
        f = REAL(-1)*cosine_function;
      if(order % 4 == 3)
        f = sine_function;
      return true;
      }
    default:
      return simplify_derivative_gen(f);
    }
  }

  template<class R, class... Args>
  bool simplify_inversion(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<INVERSION<R,Args...>>(f);
    return simplify_step(node->node);
  }
  
  template<class R, class... Args>
  bool simplify_transposition(std::shared_ptr<Node<R, Args...>>& f)
  {
    auto node = std::dynamic_pointer_cast<TRANSPOSITION<R,Args...>>(f);
    switch(node->node->get_type()){
    case ANALYTIC_OPERATION::ADDITION:
      {
        auto child = std::dynamic_pointer_cast<ADDITION<R, Args...>>(node->node);
        f = transpose(child->lhs, node->center)+transpose(child->rhs, node->center);
        return true;
      }
    case ANALYTIC_OPERATION::SUBTRACTION:
      {
        auto child = std::dynamic_pointer_cast<SUBTRACTION<R, Args...>>(node->node);
        f = transpose(child->lhs, node->center)-transpose(child->rhs, node->center);
        return true;
      }
    case ANALYTIC_OPERATION::SCALAR_ADDITION:
      {
        auto child = std::dynamic_pointer_cast<SCALAR_ADDITION<R, Args...>>(node->node);
        f = child->scalar+transpose(child->node, node->center);
        return true;
      }
    case ANALYTIC_OPERATION::MULTIPLICATION:
      {
        auto child = std::dynamic_pointer_cast<MULTIPLICATION<R, Args...>>(node->node);
        f = transpose(child->lhs, node->center)*transpose(child->rhs, node->center);
        return true;
      }
    case ANALYTIC_OPERATION::SCALAR_MULTIPLICATION:
      {
        auto child = std::dynamic_pointer_cast<SCALAR_MULTIPLICATION<R, Args...>>(node->node);
        f = child->scalar*transpose(child->node, node->center);
        return true;
      }
    case ANALYTIC_OPERATION::COMPOSITION:
      {
        auto child = std::dynamic_pointer_cast<COMPOSITION<R, Args...>>(node->node);
        f = compose(child->lhs, transpose(child->rhs, node->center));
        return true;
      }
    case ANALYTIC_OPERATION::INVERSION:
      {
        auto child = std::dynamic_pointer_cast<INVERSION<R, Args...>>(node->node);
        f = invert(transpose(child->node, node->center));
        return true;
      }
    case ANALYTIC_OPERATION::POLYNOMIAL:
      {
        f = std::dynamic_pointer_cast<poly_impl::POLY<R,Args...>>(node->node)->continuation(node->center);
        return true;
      }
    default:
      return simplify_step(node->node);
    }
  }
  
  template<class R, class... Args>
  bool simplify_step(std::shared_ptr<Node<R, Args...>>& f)
  {
    if(f->visited){
      return false;
    }
    f->visited = true;
    switch(f->get_type()){
    case ANALYTIC_OPERATION::ADDITION:
        return simplify_addition(f);
    case ANALYTIC_OPERATION::SUBTRACTION:
        return simplify_subtraction(f);
    case ANALYTIC_OPERATION::SCALAR_ADDITION:
        return simplify_scalar_addition(f);
    case ANALYTIC_OPERATION::MULTIPLICATION:
      return simplify_multiplication(f);
    case ANALYTIC_OPERATION::SCALAR_MULTIPLICATION:
      return simplify_scalar_multiplication(f);
    case ANALYTIC_OPERATION::COMPOSITION:
      return simplify_composition(f);
    case ANALYTIC_OPERATION::TRANSPOSITION:
      return simplify_transposition(f);
    case ANALYTIC_OPERATION::INVERSION:
      return simplify_inversion(f);
    case ANALYTIC_OPERATION::DERIVATIVE:
      return simplify_derivative(f);
    default:
      return false;
    }
  }

  template<class R, class... Args>
  void simplify(std::shared_ptr<Node<R, Args...>>& f)
  {
    int counter = 0;
    while(simplify_step(f)){
      ++counter;
      f->reset_visited();
    }
    f->reset_visited();
  }

  template<class R, class... Args>
  std::shared_ptr<Node<R, Args...>> prune(const std::shared_ptr<Node<R,Args...>>& node)
  {
    return std::make_shared<PRUNE<R,Args...>>(node);
  }

  template<class R, class... Args>
  int get_node_number(const std::shared_ptr<Node<R, Args...>>& node)
  {
  }
  
}
#endif

