#ifndef HYDROLOGY_H
#define HYDROLOGY_H

#include "point.hpp"

class Primitive
{
  private:
  Point loc;
  float elevation;
  int priority;
  bool isMouthNode;
  int contourIndex;

  public:
  Primitive();
  Primitive(Point loc, float elevation, int priority, bool isMouthNode, int contourIndex);
  //a trivial implicitly-declared destructor will be sufficient
};

#endif