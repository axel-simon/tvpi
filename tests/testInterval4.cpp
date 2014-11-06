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
  yBound.updateLower(-3);
  Inequality e1 = Inequality(Point(-2,0), Point(6,-3));
  Inequality e2 = Inequality(Point(6,-3), Point(-1,0));
  xBound.updateUpper(5);
  yBound.updateLower(-2);
  if (chat) cerr << "tightening upper bound : x=" << xBound << " y=" << yBound << endl;
  refineBoundsToZ<true>(xBound, yBound, &e1, &e2);
  if (chat) cerr << "after integral tightening : x=" << xBound << " y=" << yBound << endl;
  bool correct= !xBound.updateUpper(1);
  // Return 1 on error.
  return !correct;
};
