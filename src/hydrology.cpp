#include "hydrology.hpp"

Primitive::Primitive()
: loc(Point(0,0)), priority(0), isMouthNode(false), contourIndex(0)
{

}

Primitive::Primitive(Point loc, int priority, int contourIndex)
: loc(loc), priority(priority), contourIndex(contourIndex)
{

}