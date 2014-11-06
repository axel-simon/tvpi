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
  // This are two inequalities with on empty quadrant in between them. The first
  // feasible integral point lies in the upper rigth corner at (4,2).
  Inequality e1 = Inequality(Point(1,0), Point(mpq_class(7*5+3,5),4));
  Inequality e2 = Inequality(Point(mpq_class(7*5+3,5),4), Point(0,1));
  xBound.updateUpper(extremeX<true>(true, e1, e2));
  yBound.updateUpper(extremeY<true>(true, e1, e2));
  if (chat) cerr << "tightening upper bound : x=" << xBound << ", y=" << yBound  << endl;
  refineBoundsToZ<true>(xBound, yBound, &e1, &e2);
  if (chat) cerr << "after integral tightening : x=" << xBound << ", y=" << yBound << endl;
  bool correct= !xBound.updateUpper(4) && !yBound.updateUpper(2);
  // Return 1 on error.
  return !correct;
};

