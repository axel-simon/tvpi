/* Test program for removing a variable such that two adjacent terms
   can be substituted into three other equalities.
 */

#include <iostream>
#include "affine.hh"
#include "common.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Domain dom;
  DomVar vars[14];

  vector<LinComponent> emptyVec;
  Interval<true> c8 = Interval<true>();
  c8.updateUpper(4);
  c8.updateLower(-2);

  bool correct = true;
     
  // Return 1 on error.
  return !correct;
};
