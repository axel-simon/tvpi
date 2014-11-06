/* Test program for more involved interval operations.
 */

#include <iostream>
#include "interval.hh"
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Interval<true> xBound;
  Interval<true> yBound;
  Inequality e1 = Inequality(Point(-1,1), Point(6,2));
  Inequality e2 = Inequality(Point(4,1), Point(-1,3));
  xBound.updateUpper(extremeX<true>(true, e1, e2));
  if (chat) cerr << "tightening upper bound : " << xBound << endl;
  refineBoundsToZ<true>(xBound, yBound, &e1, &e2);
  if (chat) cerr << "after integral tightening : " << xBound << endl;
  // The upper x bound should be one now.
  bool correct= !xBound.updateUpper(1);
  // Return 1 on error.
  return !correct;
};
