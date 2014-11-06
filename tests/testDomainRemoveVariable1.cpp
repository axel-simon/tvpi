/* Test program for removing a variable such that two adjacent terms
   can be substituted into three other equalities.
 */

#define IN_TVPI_CHECK

#include <iostream>
#include "affine.hh"
#include "common.hh"

using namespace std;
using namespace Tvpi;

// I should fix this, this isn't testing anything.
int main(int argc) {
  Domain dom;
  DomVar vars[14];

  vector<LinComponent> emptyVec;
  Interval<true> c8 = Interval<true>();
  c8.updateUpper(4);
  c8.updateLower(-2);

  bool correct = false;
     
  // Return 1 on error.
  return !correct;
};
