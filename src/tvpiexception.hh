// tvpiexception.h
// Definition of exceptions that the TVPI library may throw.

#ifndef __TVPIEXCEPTION_H
#define __TVPIEXCEPTION_H

#include<iostream>

namespace Tvpi {
  class Exception;
  class IllegalArgument;
}; // namespace Tvpi

class Tvpi::Exception {};

//  class Unsatisfiable : TVPIException {
//   public:
//    Unsatisfiable() {};

//    friend ostream& operator<<(ostream& stream, const Unsatisfiable& e) {
//      return stream << "The polyhedron became unsatisfiable.";
//    };
//  };

class Tvpi::IllegalArgument : Tvpi::Exception {
  char* location;
  char* reason;
public:
  IllegalArgument(char* loc, char* rea) : location(loc), reason(rea) {};
  
  friend std::ostream& operator<<(std::ostream& stream,
				  const IllegalArgument& e) {
    return stream << e.location << ": Illegal Argument: " << e.reason;
  };
};
  


#endif // __TVPIEXCEPTION_H
