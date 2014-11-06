/* Test program for calculating directions of inequalities.
 */

#include <iostream>
#include "planar.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Inequality eqs[8] = {
    Inequality( 1, 0, 7),
    Inequality( 1, 1, 7),
    Inequality( 0, 1, 7),
    Inequality(-1, 1, 7),
    Inequality(-1, 0, 7),
    Inequality(-1,-1, 7),
    Inequality( 0,-1, 7),
    Inequality( 1,-1, 7)
  };
  bool correct=true;
  for (int i=0; i<8; i++) for (int j=0; j<8; j++) {

    if (chat)
      cout << "angle between " << eqs[i] << " and " << eqs[j] << " is "
	   << (eqs[i].lessThanPi(eqs[j]) ? "less than" : "greater equal")
	   << " pi" << endl;
    if (eqs[i].lessThanPi(eqs[j])!=(j-i+7)% 8 < 3) correct=false;
  };

  // Return 1 on error.
  return !correct;
};
