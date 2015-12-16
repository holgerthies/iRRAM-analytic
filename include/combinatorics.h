/*------------------------------------------------------------
 * collection of combinatorical functions
-------------------------------------------------------------*/

#ifndef COMBINATORICS_h
#define COMBINATORICS_h
#include <vector>
namespace iRRAM{
  std::vector<std::vector<unsigned long>> partitions(const unsigned long n, const unsigned long k);
  std::vector<std::vector<unsigned long>> bounded_count(const std::vector<unsigned long>& bound);
  INTEGER factorial(int n);
  INTEGER choose(int n, int k);
  REAL inv_factorial(const int n);
}
#endif
