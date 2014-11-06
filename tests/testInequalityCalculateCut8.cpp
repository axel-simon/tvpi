/* Test program for shrinking around integral grid.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  // An example the seems to loop.
  Inequality e1 = Inequality(mpz_class(-1), mpz_class(0),
			     mpz_class("2087832234"));
  Inequality e2 = Inequality(mpz_class(-3), mpz_class("-178954240"),
			     mpz_class("894672896"));
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
  bool correct=cut2->calculateCut(e2)==0;

  // Return 1 on error.
  return !correct;
};
