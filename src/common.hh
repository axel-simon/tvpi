// common.h - Common type declarations to be used by the library and
// the C wrapper.

#ifndef COMMON_H_
#define COMMON_H_

#ifdef __cplusplus
#include<cstddef>
extern "C" {
#else
#include<stddef.h>
#endif // __cplusplus

#include <limits.h>

// TVPI inequalities are never parallel to the axes and hence can take
// on four principal directions. The directions here refer to the
// normal vector of an inequality in the euclidian space. The
// directions refer to the starting angle of the normal vector. For
// instance, an inequality with positive x and y coordinate has a
// normal vector that points towards the east.
enum Direction {
  east,
  north,
  west,
  south,
  total,
  numDirs,
};

// The result of intersecting the domain with a set of inequalities.
enum Result {
  // This value is returned if the underlying domain is now empty.
  resUnsatisfiable=0,
  // This value denotes that at least some of the added inequalities
  // had an effect on the underlying domain.
  resChanged,
  // This value denotes that all passed-in inequalities had no effect
  // on the underlying domain.
  resRedundant
};

// A generic variable used in LinComponent.
typedef signed int Variable;

// An identifier number for variables in a TVPI system.
typedef size_t TvpiVar;

// An identifier number for variables in the reduced product domain.
typedef Variable DomVar;

// An abstract number denoting the type of a variable.
typedef size_t TypeId;

// The number of LSBs that are zero.
typedef unsigned char Mult;

#ifdef __cplusplus
}
#endif // __cplusplus

#endif // COMMON_H_

