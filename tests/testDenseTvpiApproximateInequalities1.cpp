/* Test program for the tvpi domain.
 *
 * Approximate an inequality with one unbounded variable.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"
#include "tvpi.hh"

using namespace std;
using namespace Tvpi;

int main(int argc) {
  bool chat = argc>1;
  DenseTvpi<true> dom(6);
  dom.createVariable(6);
  dom.createVariable(4);
  dom.createVariable(3);
  DomVar v3=dom.createVariable(2);
  DomVar v4=dom.createVariable();
  // Set v4 \in [4..6].
  Inequality* es[]= {
    new Inequality(1,1,8),
    new Inequality(-1,1,-2)
  };
  if (chat) cout << "before inserting 4 <= e <= 6\n" << dom;
  assert(dom.addInequalities(v4, v3, 2, es));
  dom.createVariable();
  if (chat) cout << "before adding inequality: " << endl << dom;
  vector<LinComponent> lc;
  //  lc.push_back(LinComponent(-1,0));
  lc.push_back(LinComponent(-2,2));
  // lc.push_back(LinComponent(3,3));
  lc.push_back(LinComponent(2,4));
  lc.push_back(LinComponent(1,5));
  int res=dom.approximateInequality(lc, -10, true);
  if (chat) cout << "after adding inequality: " << endl << dom;
  if (chat) cout << "result: " << res << endl;

  return !(res==resChanged &&
	   dom.getInterval(5).getUpper()==mpz_class(8) &&
	   dom.getInterval(5).getLower()==mpz_class(4));
};
