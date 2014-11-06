/* Test program for shrinking around integral grid.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Inequality e1 = Inequality(1, -5, 35);
  Inequality e2 = Inequality(2, -5, 39);
  assert(e1.lessThanPi(e2));
  if (chat) 
    cerr << "calculating cuts between e1=" << e1 << " and e2=" << e2 << endl;
  Inequality* cut1 = e1.calculateCut(e2);
  assert(cut1);
  if (chat) 
    cerr << "first cut is " << *cut1 << endl;
  Inequality* cut2 = cut1->calculateCut(e2);
  assert(cut2);
  if (chat) 
    cerr << "first cut is " << *cut2 << endl;
  bool correct=
    e1.lessThanPi(*cut1) &&
    cut1->lessThanPi(*cut2) &&
    cut2->lessThanPi(e2);

  // Return 1 on error.
  return !correct;
};
