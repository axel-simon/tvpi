/* Test program for shrinking around integral grid.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  // Harvey's example.
  Inequality e1 = Inequality(5, 2, 8);
  Inequality e2 = Inequality(-2, 3, 4);
  if (chat)
    cerr << "calculating cuts between e1=" << e1 << " and e2=" << e2 << endl;
  Inequality* cut1 = e1.calculateCut(e2);
  assert(cut1);
  if (chat)
    cerr << "first cut is " << *cut1 << endl;
  Inequality* cut2 = cut1->calculateCut(e2);
  assert(cut2);
  if (chat)
    cerr << "second cut is " << *cut2 << endl;
  Inequality* cut3 = cut2->calculateCut(e2);
  assert(cut3);
  if (chat)
    cerr << "third cut is " << *cut3 << endl;
  bool correct=cut3->calculateCut(e2)==0;

  // Return 1 on error.
  return !correct;
};
