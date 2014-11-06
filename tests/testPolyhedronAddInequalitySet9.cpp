/* Bug in redundancy removal where two inequalities that are more than pi
   appart do not tighten the bounds completely.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Polyhedron<true> p;
 
  vector<Inequality*> eqs;
  eqs.push_back(new Inequality(-4079,4080,4080));
  eqs.push_back(new Inequality(4079,-4080,-4080));
  Interval<true> xBound;
  Interval<true> yBound;
  xBound.updateLower(2041);
  xBound.updateUpper(4080);
  yBound.updateLower(2041);
  yBound.updateUpper(4080);
  p.addInequalitySet(xBound, yBound, eqs);
  if (chat) cout << "result:\n" << p
		 << "x in " << xBound << ", y in " << yBound << endl;

  // The inequalities only have a single integral point as solution.
  bool correct=!xBound.updateLower(4080) && !xBound.updateUpper(4080);

  // Return 1 on error.
  return !correct;
};
