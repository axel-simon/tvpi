// The the continued fraction functions.

#include<iostream>
#include<planar.hh>

using namespace std;
using namespace Tvpi;

int main(int argc) {
  bool chat = argc>1;

  // Example fractions from the Davenport book.
  int a[] = { 2, 3, 67, 24, 42, 31 };
  int b[] = { 3, 2, 24, 67, 31, 42 };
	 
  bool correct = true;
  Convergents c;

  for (size_t i=0; i<sizeof(a)/sizeof(a[0]); i++) {
    mpz_class aa(a[i]);
    mpz_class bb(b[i]);
  
    size_t last = c.initialize(aa,bb)-1;
    if (chat) cerr << c << endl;

    // Get the next-to-last approximation.
    const Convergents::Approximation* app = c[last-1];

    if (last & 1) {
      // last is odd, this equation is 1
      if (aa*app->B - bb*app->A !=  1) correct=false;
      if (chat)	cerr << aa << " * " << app->B << " - "
		     << bb << " * " << app->A << " == "
		     << aa*app->B - bb*app->A << endl;
    } else {
      // last is even, this equation is -1
      if (aa*app->B - bb*app->A != -1) correct=false;
      if (chat)	cerr << aa << " * " << app->B << " - "
		     << bb << " * " << app->A << " == "
		     << aa*app->B - bb*app->A << endl;
    };
  };

  // Return 1 on error.
  return !correct;
};
