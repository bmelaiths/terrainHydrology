#ifndef HYDROLOGY_H
#define HYDROLOGY_H

#include <vector>
#include <stdio.h>
#include <stdint.h>

#include "point.hpp"
#include "kdtree.hpp"

class Primitive
{
  public:
  size_t id, parent;
  std::vector<size_t> children;
  bool isMouthNode;
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
  size_t binarySize();
  void toBinary(uint8_t *buffer);
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
  KDTree<size_t> tree;

  public:
  Primitive addMouthNode(
    Point loc, float elevation, int priority, int contourIndex
  );
  Primitive addRegularNode(
    Point loc, float elevation, int priority, size_t parent
  );
  std::vector<Edge> edgesWithinRadius(Point loc, float radius);
  Primitive getNode(size_t idx);
};

void writeBinary(Hydrology hydrology, FILE *stream);

#endif