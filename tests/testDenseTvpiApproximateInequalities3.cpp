/* Test program for the tvpi domain.
 *
 * Approximate an inequality with no unbounded variable.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"
#include "tvpi.hh"

using namespace std;
using namespace Tvpi;

int main(int argc) {
  bool chat = argc>1;
  TvpiVar ident[] = {0, 1, 2, 3, 4, 5, 6};
  DenseTvpi<false> dom(6);
  DomVar v1=dom.createVariable(mpz_class(30), 1);
  DomVar v2=dom.createVariable(mpz_class(15), 2);
  DomVar v3=dom.createVariable(mpz_class(10), 3);
  DomVar v4=dom.createVariable(Interval<false>(),4);
  // Set v4 \in [20..30].
  Inequality* es[]= {
    new Inequality(1,1,40),
    new Inequality(1,-1,-10)
  };
  assert(dom.addInequalities(v3, v4, 2, es));
  DomVar v5=dom.createVariable(Interval<false>(),5);
  if (chat) cerr << "before doing anything:" << endl << dom;
  vector<LinComponent> lc;
  lc.push_back(LinComponent(-1,v1));
  lc.push_back(LinComponent(-2,v2));
  lc.push_back(LinComponent(3,v3));
  lc.push_back(LinComponent(2,v4));
  lc.push_back(LinComponent(-1,v5));
  DenseTvpi<false> dom2=dom;
  bool correct = dom.approximateInequality(lc, 55, true);
  // The value of of the linear expression is between 65 and 85,
  // depending if we minimize or maximize the expression.
  if (chat) cerr << "after inserting big inequality:" << endl << dom;
  correct= dom2.includes(dom, ident) && (!dom.includes(dom2, ident))
    && correct;
  // Inverted logic for the return value.
  return !correct;
};
