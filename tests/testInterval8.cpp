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
  // This is an equality.
  Inequality e1 = Inequality(Point(1,0), Point(mpq_class(7*5+3,5),4));
  Inequality e2 = Inequality(Point(mpq_class(7*5+3,5),4), Point(1,0));
  xBound.updateUpper(7);
  yBound.updateUpper(4);
  xBound.updateLower(-100);
  yBound.updateLower(-200);
  if (chat) cerr << "tightening upper bound : x=" << xBound << ", y=" << yBound  << endl;
  refineBoundsToZ<true>(xBound, yBound, &e1, &e2);
  if (chat) cerr << "after tightening upper: x=" << xBound << ", y=" << yBound << endl;
  refineBoundsToZ<true>(xBound, yBound, &e2, &e1);
  if (chat) cerr << "after tigthening lower : x=" << xBound << ", y=" << yBound << endl;
  bool correct= !xBound.updateUpper(1) && !yBound.updateUpper(0) &&
		!xBound.updateLower(-98) && !yBound.updateLower(-60);
  // Return 1 on error.
  return !correct;
};

