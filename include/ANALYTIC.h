/*-------------------------------------------------
 * Class for real analytic functions on [0,1]
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
#include <vector>
namespace iRRAM
{
  class ANALYTIC{
    private:
      int l;
      std::vector<std::shared_ptr<POWERSERIES<1,REAL>>> pwrs;
    public:
      ANALYTIC(std::shared_ptr<POWERSERIES<1,REAL>>, const int l);
      ANALYTIC(std::shared_ptr<std::function<REAL(const std::vector<unsigned long>&)>> f, const int l):
        ANALYTIC(std::shared_ptr<POWERSERIES<1,REAL>>(new POWERSERIES<1,REAL>(f)), l) {};
      REAL eval_k(const REAL& x, const int k);
      REAL operator ()(const REAL&) const;
      int get_known_coeffs(const int k) const;
      int get_l() const;
      std::shared_ptr<POWERSERIES<1,REAL>> get_series(const int k) const;
  };
  // compute the d-th derivative of the analytic function f
  ANALYTIC derive(const ANALYTIC& f, const int d);
}
#endif
