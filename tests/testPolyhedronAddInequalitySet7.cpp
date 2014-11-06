/* Test program for tightening.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Polyhedron<true> p;
 
  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(-34,55,1));
  Interval<true> xBound;
  xBound.updateUpper(0);
  Interval<true> yBound;
  p.addInequalitySet(xBound, yBound, eqs);
  //cout << "before:\n" << p
  //     << "x in " << xBound << ", y in " << yBound << endl;

  bool correct=p.getNoOfInequalities()==2;

  // Return 1 on error.
  return !correct;
};
