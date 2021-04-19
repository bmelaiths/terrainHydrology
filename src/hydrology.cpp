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

Edge::Edge(Primitive node0, Primitive node1)
: node0(node0), node1(node1)
{
}

bool ComparePrimitive::operator() (const Primitive& a, const Primitive& b) const {
  return a.priority < b.priority;
}