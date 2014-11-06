/* Test program for basic interval operations.
 */

#include <iostream>
#include "interval.hh"

using namespace std;
using namespace Tvpi;

  
int main() {
  Interval<true> i1;
  i1.updateUpper(8);
  Interval<true> i2;
  i2.updateLower(4);
  Interval<true> i3=i2;
  i3.updateUpper(8);
  Interval<true> i4;
  bool correct= i1.includes(i3) && i2.includes(i3) &&
    i4.includes(i1) && i4.includes(i2) && i4.includes(i3);

  // Return 1 on error.
  return !correct;
};
