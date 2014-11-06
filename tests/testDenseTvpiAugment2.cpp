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

  DenseTvpi<true> dom;
  TvpiVar x = dom.createVariable(0);
  TvpiVar y = dom.createVariable(1);

  dom.augment(y,x);
 
  if (chat) cout << "result after augmenting x=0 with y=1:" << endl << dom;

  bool correct =
    dom.minValue(x)==0 && dom.maxValue(x)==1 &&
    dom.minValue(y)==1 && dom.maxValue(y)==1;

  return !correct;
};
