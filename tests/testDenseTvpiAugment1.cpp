/* Test program for additive update.
 */

#include <iostream>
#include "planar.hh"
#include "polyhedron.hh"
#include "tvpi.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;

  TvpiVar ident[] = {0,1,2,3,4};
  const TvpiVar maxVars=5;
  DenseTvpi<true> dom1;
  TvpiVar vars[maxVars];
  vars[0]=dom1.createVariable();
  for (TvpiVar i=1; i<maxVars; i++) {
    vars[i]=dom1.createVariable();
    vector<Inequality*> eqs;
    eqs.push_back(new Inequality(-1, 1, 1));
    eqs.push_back(new Inequality(1, -1, -1));
    dom1.addInequalitySet(vars[i-1], vars[i], eqs);
  };
  if (chat) cout << "initial: " << endl;
  if (chat) dom1.showDist(cout) << dom1;
  DenseTvpi<true> dom2=dom1;
  dom2.augment(2,4);
  if (chat) cout << "before:" << endl;
  if (chat) dom1.showDist(cout) << dom1;
  if (chat) cout << "augmenting 4 with 2:" << endl;
  if (chat) dom2.showDist(cout) << dom2;
  if (chat) cout << "dom1.includes(dom2) " << dom1.includes(dom2, ident)
		 << endl;
  if (chat) cout << "dom2.includes(dom1) " << dom2.includes(dom1, ident)
		 << endl;
  bool correct=
    !dom1.includes(dom2, ident) &&
    dom2.includes(dom1, ident);
  return !correct;
};
