/*-------------------------------------------------
 * Class for real analytic functions on [0,1]
 ------------------------------------------------*/
#ifndef ANALYTIC_H
#define ANALYTIC_H
#include "POWERSERIES.h"
class ANALYTIC{
  private:
    int l;
    vector<POWERSERIES<1,REAL>> pwrs;
  public:
    ANALYTIC(const std::function<int,REAL>& f, const int l);
    REAL operator ()(const REAL&) const;
}
#endif
