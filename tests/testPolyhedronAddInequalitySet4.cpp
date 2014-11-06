/* Test program to check redundancy removal.
 */

#include <iostream>
#include "interval.hh"
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Polyhedron<true> p;

  // Create an equality.
  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(Point(13,-1), Point(11,4)));
  eqs.push_back(new Inequality(Point(11,4), Point(13,-1)));

  Interval<true> xBound;
  Interval<true> yBound;
  xBound.updateLower(mpq_class(10));
  xBound.updateUpper(mpq_class(10));
  p.addInequalitySet(xBound, yBound, eqs);
  if (chat) cout << p << "x: " << xBound << " y: " << yBound << endl;

  bool correct=p.getNoOfInequalities()==0;
  // Return 1 on error.
  return !correct;
};
