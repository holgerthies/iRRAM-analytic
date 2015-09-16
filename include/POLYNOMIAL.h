/*-------------------------------------
 * Class for multivariate polynomials
 -------------------------------------*/
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H
#include <vector>
#include <memory>
#include <string>
template<unsigned int n, class T>
class POLYNOMIAL{
  private:
    std::vector<std::shared_ptr<POLYNOMIAL<n-1,T>>> coefficients;
  public:
    // empty constructor for 0 polynomial
    POLYNOMIAL() : coefficients(std::vector<std::shared_ptr<POLYNOMIAL<n-1,T>>>{std::shared_ptr<POLYNOMIAL<n-1,T>>(new POLYNOMIAL<n-1,T>())}) {};
    
    std::string to_string(){
      char letter='w'+n;
      std::string result;
      return result;
    }
};

// Template specialization for univariate polynomials
template<class T>
class POLYNOMIAL<1,T>{
  private:
    std::vector<T> coefficients;
  public:
    // empty constructor for 0 polynomial
    POLYNOMIAL() : coefficients(std::vector<T>{T()}) {};
    // constructor with given coefficients
    POLYNOMIAL(std::vector<T> coeffs) : coefficients(coeffs) {};
    // Horner Scheme for evaluation
    T operator ()(const T& x) const{
      T ans=coefficients.back();
      for(int i=coefficients.size()-2; i>=0; i--){
        ans = coefficients[i]+ans*x;
      }
      return ans;
    }
};
#endif
