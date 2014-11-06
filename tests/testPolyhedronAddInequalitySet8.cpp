/* Bug in redundancy removal where a single or pi-apart inequality
   didn't propagate a bound.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Polyhedron p;
 
  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(-1,2,0));
  eqs.push_back(new Inequality(1,-2,5));
  Interval<false> xBound;
  Interval<false> yBound;
  xBound.updateLower(7);
  xBound.updateUpper(12);
  p.addInequalitySet(xBound, yBound, eqs);
  if (chat) cout << "result:\n" << p
		 << "x in " << xBound << ", y in " << yBound << endl;

  bool correct=!yBound.updateLower(1) && !yBound.updateUpper(6);

  // Return 1 on error.
  return !correct;
};
