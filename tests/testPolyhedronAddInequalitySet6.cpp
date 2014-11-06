/* Test program for tightening. Used in paper.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Polyhedron<true> p;
 
  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(Point(9,1), Point(2,9)));
  eqs.push_back(new Inequality(Point(13,6), Point(2,3)));
  eqs.push_back(new Inequality(Point(2,3), Point(3,-1)));
  eqs.push_back(new Inequality(Point(-2,-2), Point(9,3)));
  Interval<true> xBound;
  Interval<true> yBound;
  p.addInequalitySet(xBound, yBound, eqs);
  //cout << "before:\n" << p
  //     << "x in " << xBound << ", y in " << yBound << endl;

  bool correct=p.getNoOfInequalities()==5;

  // Return 1 on error.
  return !correct;
};
