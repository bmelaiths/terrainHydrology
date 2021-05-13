#ifndef HYDROLOGY_H
#define HYDROLOGY_H

#include <vector>
#include <list>
#include <stdio.h>
#include <stdint.h>

#include "point.hpp"
#include "kdtree.hpp"
#include "forest.hpp"

class Hydrology;

class Primitive
{
  private:
  size_t id;
  Primitive *parent;
  std::vector<Primitive*> children;
  Point loc;
  float elevation;
  int priority;
  int contourIndex;

  public:
  friend Hydrology;

  public:
  Primitive();
  Primitive
  ( //for mouth nodes
    size_t id, Point loc, float elevation, int priority, int contourIndex
  );
  Primitive
  ( //for regular nodes
    size_t id, Primitive *parent, Point loc, float elevation, int priority
  );
  //a trivial implicitly-declared destructor will be sufficient
  size_t binarySize();
  void toBinary(uint8_t *buffer);

  size_t getID();
  Primitive* getParent();
  bool hasParent();
  std::vector<Primitive*> getChildren();
  size_t numChildren();
  Point getLoc();
  float getElevation();
  int getPriority();
  int getContourIndex();
};

class Edge
{
  public:
  Primitive *node0, *node1;

  public:
  Edge(Primitive *node0, Primitive *node1);
};

class Hydrology
{
  private:
  omp_lock_t lock;

  public:
  std::vector<Primitive*> indexedNodes;
  Forest<Primitive*> trees;

  public:
  Hydrology();
  Hydrology(Point lowerLeft, Point upperRight);
  ~Hydrology();
  Hydrology(const Hydrology& other);
  Hydrology(Hydrology&& other);

  Hydrology& operator=(const Hydrology& other);
  Hydrology& operator=(Hydrology&& other);

  void lockNetwork();
  Primitive* addMouthNode(
    Point loc, float elevation, int priority, int contourIndex
  );
  Primitive* addRegularNode(
    Point loc, float elevation, int priority, size_t parent
  );
  void unlockNetwork();
  AreaLock lockArea(Point loc, float radius);
  std::vector<Edge> queryArea(Point loc, float radius);
  Primitive getNode(size_t idx);
  size_t numNodes();
};

void writeBinary(Hydrology& hydrology, FILE *stream);

#endif