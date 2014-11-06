/* Test program for basic interval operations.
 */

#include <iostream>
#include "interval.hh"

using namespace std;
using namespace Tvpi;

  
int main(int argc) {
  bool chat = argc>1;
  Interval<true> i1=Interval<true>(7);
  if (chat) cout << "i1: " << i1 << endl;
  Interval<true> i2=Interval<true>(16);
  if (chat) cout << "i2: " << i2 << endl;
  Interval<true> i3=Interval<true>(i1,i2);
  if (chat) cout << "i3: " << i3 << endl;
  Interval<true> i4=i1;
  i4.widen(i3, mpz_class(0));
  if (chat) cout << "i1 widened by i3: " << i4 << endl;
  i2.updateUpper(8);
  if (chat) cout << "i2: " << i2 << " is empty: " << i2.isEmpty() << endl;
  if (chat) cout << "i1 |= i3: " << i3.includes(i1) << endl;
  bool correct= i2.isEmpty() && i3.includes(i1);

  // Return 1 on error.
  return !correct;
};
