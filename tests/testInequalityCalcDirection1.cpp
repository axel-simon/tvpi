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
  for (int i=0; i<8; i++) {
    if (chat)
      cout << "dir of " << eqs[i] << " : " << eqs[i].calcDirection() << endl;
    if (eqs[i].calcDirection()!=(i/2)) correct=false;
  };

  // Return 1 on error.
  return !correct;
};
