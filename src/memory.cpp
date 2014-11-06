#include "memory.hh"

#ifdef DEBUG_MEMORY

#include <iostream>

namespace Tvpi {
  std::ofstream o("allocations.txt");
  size_t allocIdx = 0;
}


#endif
