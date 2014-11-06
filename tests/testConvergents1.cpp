// The the continued fraction functions.

#include<iostream>
#include<planar.hh>

using namespace std;
using namespace Tvpi;

int main(int argc) {
  bool chat = argc>1;

  Convergents c;
  c.initialize(12,19);
  
  for (size_t i = 18; i>0; i--) {
    mpz_class i_(i);
    if (chat)
      cerr << "convergent smaller than " << i << ": "
	   << c.seekClosest(i_)
	   << endl << c << endl;
  };

  c.initialize(12,19);

  mpz_class one(1);
  mpz_class two(2);
  mpz_class five(5);
  mpz_class eight(8);
  bool correct =
    c.seekClosest(eight)==mpq_class(5,8) &&
    c.seekClosest(five)==mpq_class(3,5) &&
    c.seekClosest(two)==mpq_class(1,2) &&
    c.seekClosest(one)==0;

  // Return 1 on error.
  return !correct;
};
