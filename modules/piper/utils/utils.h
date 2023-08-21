#ifndef PIPER_UTILS_H
#define PIPER_UTILS_H

namespace piper
{

enum class ComponentCategory : int
{
  BoundaryLike = 0,
  Volumetric = 1,
  JunctionLike = 2
};

}

#endif // PIPER_UTILS_H
