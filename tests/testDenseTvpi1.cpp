/* Test program for approximateInequalities
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"
#include "tvpi.hh"

using namespace std;
using namespace Tvpi;

// Insert inequalities that say that var1 and var2 are equal.
void makeEqual(DenseTvpi& tvpi, Variable var1, Variable var2) {
  Inequality* eq[4] = {
    new Inequality(1, -1, 0),
    new Inequality(-1, 1, 0)
  };
  tvpi.addInequalitySet(var1,
			var2,
			2,
			&eq[0],
			false);
}
  
int main() {
  const int maxVars=10;
  DenseTvpi dom(maxVars);
  size_t vars[maxVars];
  for (int i=0; i<maxVars/2; i++) vars[i]=dom.createVariable(i);
  for (int i=maxVars/2+1; i<maxVars; i++) vars[i]=dom.createVariable();
  cout << "before projection: " << endl << dom;
  size_t proj[] = { 0, 6, 4, 7, 3, 5 };
  dom.projectOnto(sizeof(proj)/sizeof(proj[0]), proj);
  cout << "after projection: " << endl << dom;
};
