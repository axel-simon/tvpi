/* Test program for shrinking around integral grid.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Inequality e1 = Inequality(3,1,25);
  Inequality e2 = Inequality(3,5,46);
  if (chat)
    cerr << "calculating cuts between e1=" << e1 << " and e2=" << e2 << endl;
  Inequality* cut1 = e1.calculateCut(e2);
  assert(cut1);
  Inequality* cut2 = cut1->calculateCut(e2);
  assert(cut2);
  Inequality* cut3 = cut2->calculateCut(e2);
  assert(!cut3);
  bool correct=
    !cut1->includes(e1,*cut2) &&
    !cut2->includes(*cut1,e2);

  // Return 1 on error.
  return !correct;
};
