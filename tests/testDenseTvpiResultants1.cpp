/* Test program for the closure and resultants calculation.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"
#include "tvpi.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  DenseTvpi<false> dom;
  DomVar x = dom.createVariable(1);
  DomVar y = dom.createVariable(2);
  DomVar z = dom.createVariable();
  Inequality* eqsXZ[] = {
    new Inequality(2, -3, 5),
    new Inequality(-2, 4, -5)
  };
  dom.addInequalities(x, z, 2, eqsXZ);
  //cout << dom;
  // This renders the polyhedron unsatisfiable, but it doesn't really
  // exercise the closure code.
  Inequality* eqsXY[] = {
    new Inequality(1, -1, -2)
  };

  // That inequality should render the system unfeasible.
  bool correct = !dom.addInequalities(x, y, 1, eqsXY);

  return !correct;

};
