/* Test program for the tvpi domain.
 *
 * Approximate an inequality with two unbounded variable.
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
  DomVar x = dom.createVariable();
  Interval<true> greaterSix;
  greaterSix.updateUpper(6);
  Interval<true> blah;
  bool sat = dom.intersectBound(x, greaterSix);
  assert(sat);
  dom.createVariable(4);
  DomVar y = dom.createVariable();
  Interval<true> greaterThree;
  greaterThree.updateUpper(3);
  sat = dom.intersectBound(y, greaterThree);
  assert(sat);
  dom.createVariable(2);
  dom.createVariable();
  dom.createVariable();
  if (chat) cout << "before adding inequality: " << endl << dom;
  vector<LinComponent> lc;
  lc.push_back(LinComponent(-1,0));
  lc.push_back(LinComponent(-2,2));
  lc.push_back(LinComponent(2,4));
  lc.push_back(LinComponent(1,5));
  int res=dom.approximateInequality(lc, 10, false);
  if (chat) cout << "after adding inequality: " << endl << dom;
  if (chat) cout << "result: " << res << endl;
  return !(res==1);
};
