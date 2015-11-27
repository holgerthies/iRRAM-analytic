#include "iRRAM.h"
#include "POWERSERIES.h"
#include <iostream>
using namespace std;
using namespace iRRAM;
void compute(){
  using seq = std::function<POWERSERIES<1,REAL>(const unsigned long)>;
  using seq_ptr = shared_ptr<seq>;
  seq series = [] (const unsigned long n) {
    auto subseq = [n] (unsigned long m) {return REAL(int(n+m));};
    return POWERSERIES<1,REAL>(subseq);
	  
  };

  function<REAL(unsigned long, unsigned long, unsigned long)> series2 = [] (unsigned long n, unsigned long m, const unsigned long k) {return REAL(int(n+m-k));};

  shared_ptr<POWERSERIES<2,REAL>> test(new POWERSERIES<2,REAL>(seq_ptr(new seq(series))));
  POWERSERIES<3,REAL> test2(series2);
  //iRRAM::cout << evaluate_partial(test, {3,5}, 2, 5) << endl;
  iRRAM::cout << (*test)[300][245] << endl;
  iRRAM::cout << test2[300][245][1000] << endl;
  iRRAM::cout << ((*test)*(*test))[1][1] << endl;
}

