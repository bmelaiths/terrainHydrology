#include "hydrology.hpp"

Primitive::Primitive()
:
loc(Point(0,0)), elevation(0.0f), priority(0),
isMouthNode(false), contourIndex(0)
{

}

Primitive::Primitive
(
  Point loc, float elevation, int priority,
  bool isMouthNode, int contourIndex
)
:
loc(loc), elevation(elevation), priority(priority),
isMouthNode(isMouthNode), contourIndex(contourIndex)
{

}