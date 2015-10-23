/*------------------------------------------------------------
 * collection of combinatorical functions
-------------------------------------------------------------*/

#ifndef COMBINATORICS_h
#define COMBINATORICS_h
#include <vector>
namespace iRRAM{
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k);
  INTEGER factorial(int n);
  REAL choose(int n, int k);
}
#endif
