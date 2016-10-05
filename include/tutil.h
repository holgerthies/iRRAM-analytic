#ifndef TUTIL_H
#define TUTIL_H
#include <utility>
  
namespace tutil{

  template <typename T>
  T get(int index)
  {
    return T();
  }

  template <typename T, typename... Ts>
  T get(int index, T first, Ts... rest)
  {
    if(index == 0){
      return first;
    }
    return get<T>(index-1, rest...);
  }

  template <size_t n, template<typename...> class U, typename T, typename... Ts>
  struct repeat_impl
  {
    using type = typename repeat_impl<n-1,U, T, T, Ts...>::type;
  };

  template <template<typename...> class U, typename T, typename... Ts>
  struct repeat_impl<1, U, T, Ts...>
  {
    using type = U<T, Ts...>;
  };

  template <size_t n, template<typename...> class U, typename T>
  using repeat = typename repeat_impl<n, U, T>::type;
  
  template <size_t n, typename T>
  using n_tuple = repeat<n, std::tuple, T>;

  
  template <size_t n, template<typename> class U, typename T>
  struct concat_impl
  {
    using type = typename concat_impl<n-1,U, U<T>>::type;
  };

  template <template<typename> class U, typename T>
  struct concat_impl<0, U, T>
  {
    using type = T;
  };

  template <size_t n, template<typename> class U, typename T>
  using concat = typename concat_impl<n, U, T>::type;

  template <typename T, typename... Ts>
  struct sum_tuple_helper
  {
    template<typename... Us>
    static std::tuple<T, Ts...> get(const std::tuple<Us...>& lhs, const std::tuple<Us...>& rhs)
    {
      return std::tuple_cat(sum_tuple_helper<Ts...>::get(lhs, rhs), std::make_tuple(std::get<sizeof...(Ts)>(lhs)+std::get<sizeof...(Ts)>(rhs)));
    }
  };

  template <typename T>
  struct sum_tuple_helper<T>
  {
    template<typename... Us>
    static std::tuple<T> get(const std::tuple<Us...>& lhs, const std::tuple<Us...>& rhs)
    {
      return std::make_tuple<T>(std::get<0>(lhs)+std::get<0>(rhs));
    }
  };

  template <size_t... Is>
  struct seq {};

  template <size_t N, size_t M, size_t... Is>
  struct gen_seq : gen_seq<N, M-1, M-1, Is...>{};

  template <size_t N, size_t... Is>
  struct gen_seq<N, N, Is...> : seq<Is...> {};

  template <typename... Ts>
  std::tuple<Ts...> sum_tuples(const std::tuple<Ts...>& lhs, const std::tuple<Ts...>& rhs)
  {
    return sum_tuple_helper<Ts...>::get(lhs, rhs);
  }

  template <typename... Ts, size_t... Is>
  auto get_index(const std::tuple<Ts...>& tuple, seq<Is...> S) -> decltype(std::make_tuple(std::get<Is>(tuple)...)){
    return std::make_tuple(std::get<Is>(tuple)...);
  }

  template <size_t I, typename T, typename... Ts>
  std::tuple<Ts...> tuple_replace(const std::tuple<Ts...>& tuple, const T& v)
  {
    return std::tuple_cat(get_index(tuple, gen_seq<0,I>()), std::make_tuple(v), get_index(tuple, gen_seq<I+1, sizeof...(Ts)>()));
  }

  template <typename T, typename... Ts>
  std::tuple<Ts...> tail(const std::tuple<T, Ts...>& tuple){
    return get_index(tuple, gen_seq<1,sizeof...(Ts)+1>());
  }

  template<typename T>
  std::string to_string(const std::tuple<T>& tuple){
    return std::to_string(std::get<0>(tuple));
  }

  template<typename... Ts>
  std::string to_string(const std::tuple<Ts...>& tuple){
    return std::to_string(std::get<0>(tuple))+", "+to_string(tail(tuple));
  }

}
#endif
