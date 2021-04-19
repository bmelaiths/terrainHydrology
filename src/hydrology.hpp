#ifndef HYDROLOGY_H
#define HYDROLOGY_H

#include <vector>

#include "point.hpp"

class Primitive
{
  public:
  int id;
  Point loc;
  float elevation;
  int priority;
  bool isMouthNode;
  int contourIndex;
  int parent;

  public:
  Primitive();
  Primitive
  (
    Point loc, float elevation, int priority,
    bool isMouthNode, int contourIndex
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

  public:
  Hydrology(/* args */);
  ~Hydrology();

  Primitive addNode(
    Point loc, float elevation, int priority, int parent
  );
  std::vector<Edge> edgesWithinRadius();
};

class ComparePrimitive
{
  public:
  bool operator() (const Primitive& a, const Primitive& b) const;
};

#endif