/* Test program to check substitution of variables with equal variables.
 */

#include <iostream>
#include "affine.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;

  const TvpiVar maxVars=5;
  Domain dom1;
  dom1.createVariable();
  for (TvpiVar i=1; i<maxVars; i++) {
    dom1.createVariable();
    vector<LinComponent> comps;
    // Say that x_2 + x_3 = 20
    if (i==3) {
      comps.push_back(LinComponent(mpz_class(1), i-1));
      comps.push_back(LinComponent(mpz_class(1), i));
      Result res = dom1.inequality(comps, mpz_class(20), 0, true);
      assert(res==resChanged);
    } else {
      comps.push_back(LinComponent(mpz_class(1), i-1));
      comps.push_back(LinComponent(mpz_class(-2), i));
      Result res = dom1.inequality(comps, mpz_class(0), 0, true);
      assert(res==resChanged);
    }
  };
  if (chat) cout << "initial: " << endl;
  if (chat) cout << dom1;

  // Add something that has equalities in it.
  vector<LinComponent> comps;
  comps.push_back(LinComponent(mpz_class(1), 0));
  comps.push_back(LinComponent(mpz_class(1), 1));
  comps.push_back(LinComponent(mpz_class(1), 2));
  comps.push_back(LinComponent(mpz_class(1), 3));
  comps.push_back(LinComponent(mpz_class(1), 4));
  Result res = dom1.inequality(comps, mpz_class(100));

  if (chat) cout << "after intersection: " << endl;
  if (chat) cout << dom1;
   
  bool correct= res==resChanged;
  return !correct;
};
