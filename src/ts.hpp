#ifndef TERRAINS_H
#define TERRAINS_H

#include <uchar.h>
#include <vector>

#include "point.hpp"

class T
{
private:
  Point loc;
  size_t cellID;
  float elevation;
public:
  T(Point loc, size_t cellID)
  : loc(loc), cellID(cellID)
  {};
  const Point& getLoc() const {return loc;}
  size_t getCellID() const {return cellID;}
  float getElevation() const {return elevation;}
  void setElevation(float newElevation) {elevation = newElevation;}
};

class Terrain
{
private:
  std::vector<T*> allTs;
public:
  ~Terrain();
  void dumpT(Point loc, size_t cellID);
  size_t numTs() const;
  T& getT(size_t id);
};

#endif