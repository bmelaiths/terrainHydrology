#ifndef HYDROLOGY_H
#define HYDROLOGY_H

#include <vector>

#include "point.hpp"
#include "kdtree.hpp"

class Primitive
{
  public:
  size_t id, parent, child;
  bool isMouthNode, hasChild = false;
  Point loc;
  float elevation;
  int priority;
  int contourIndex;

  public:
  Primitive();
  Primitive
  ( //for mouth nodes
    size_t id, Point loc, float elevation, int priority, int contourIndex
  );
  Primitive
  ( //for regular nodes
    size_t id, size_t parentID, Point loc, float elevation, int priority
  );
  //a trivial implicitly-declared destructor will be sufficient
};

class Edge
{
  public:
  Primitive node0, node1;

  public:
  Edge(Primitive node0, Primitive node1);
};

class Hydrology
{
  private:
  /* data */

  public:
  std::vector<Primitive> indexedNodes;
  KDTree tree;

  public:
  Primitive addMouthNode(
    Point loc, float elevation, int priority, int contourIndex
  );
  Primitive addRegularNode(
    Point loc, float elevation, int priority, size_t parent
  );
  std::vector<Edge> edgesWithinRadius(Point loc, float radius);
};

class ComparePrimitive
{
  public:
  bool operator() (const Primitive& a, const Primitive& b) const;
};

#endif