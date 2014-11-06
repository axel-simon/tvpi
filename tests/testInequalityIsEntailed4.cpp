/* Test program for the planar entailment.
 */

#include <iostream>
#include "planar.hh"
#include "interval.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  // From the convex hull example.
  Inequality eq1 = Inequality(-1, 1, 10);
  Inequality et = Inequality(-4, 11, 68);
  Interval<false> xBound;
  Interval<false> yBound;
  xBound.updateUpper(mpq_class(-1));
  xBound.updateLower(mpq_class(-8));
  yBound.updateUpper(mpq_class(4));
  yBound.updateLower(mpq_class(-1));
  //  cerr << et << " includes " << eq1 << " and x in " << xBound << ", y in "
  //       << yBound << ": " << et.includes(eq1, xBound, yBound) << endl;
  bool correct=et.includes<false>(eq1, xBound, yBound);
  // Return 1 on error.
  return !correct;
};
