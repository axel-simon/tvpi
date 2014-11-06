/* Test program for shrinking around integral grid.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat=argc>1;
  mpz_class a("38326571");
  mpz_class c("156372410812");
  Inequality* cut = new Inequality(a, 1132, c);
  Inequality e = Inequality(3, 1, 12242);
  int count=0;
  while (cut) {
    count++;
    if (chat) {
      Point p(e, *cut);
      cerr << "e=" << e << " and cut=" << *cut << " intersect at "
	   << p << endl;
    };
    Inequality* newCut = cut->calculateCut(e);
    delete cut;
    cut=newCut;
  };
  if (chat) cerr << "calculated " << count << " cuts" << endl;
  bool correct=count>0;

  // Return 1 on error.
  return !correct;
};
