/* Test program for removal of inequalities when tightening bounds.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Polyhedron p;
 
  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(Point(13,-6), Point(12,-1)));
  eqs.push_back(new Inequality(Point(12,-1), Point(11,2)));
  eqs.push_back(new Inequality(Point(11,2), Point(9,5)));
  eqs.push_back(new Inequality(Point(9,5), Point(2,8)));
  eqs.push_back(new Inequality(Point(0,8), Point(-7,7)));
  eqs.push_back(new Inequality(Point(-7,7), Point(-10,6)));
  eqs.push_back(new Inequality(Point(-10,6), Point(-12,3)));
  eqs.push_back(new Inequality(Point(-12,3), Point(-13,-4)));
  eqs.push_back(new Inequality(Point(-13,-4), Point(-8,-7)));

  Interval<false> xBound;
  Interval<false> yBound;
  p.addInequalitySet(xBound, yBound, eqs);
  //  cout << "before:\n" << p
  //       << "x in " << xBound << ", y in " << yBound << endl;
  yBound.updateUpper(5);
  yBound.updateLower(-4);
  p.propagateYBounds(xBound, yBound);
  p.enforceBounds(xBound, yBound);
  //  cout << "after:\n" << p
  //       << "x in " << xBound << ", y in " << yBound << endl;
  bool correct=p.getNoOfInequalities()==5;

  // Return 1 on error.
  return !correct;
};
