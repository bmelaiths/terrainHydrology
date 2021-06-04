#ifndef HONEYCOMB_H
#define HONEYCOMB_H

#include <vector>
#include <map>

#include "point.hpp"

/**
 * @brief Represents a ridge point
 * 
 */
class Q
{
private:
  Point position;
  float elevation;
  size_t vorIndex;
  std::vector<size_t> nodes;
public:
  Q
  (
    Point position, float elevation,
    size_t vorIndex, std::vector<size_t> nodes
  ):
  position(position), elevation(elevation), vorIndex(vorIndex),
  nodes(nodes)
  {}
  ~Q() = default;
  const Point getPosition() const {return position;}
  float getElevation() const {return elevation;}
  size_t getVorIndex() const {return vorIndex;}
  std::vector<size_t> getNodes() const {return nodes;}
};

class Ridge
{
private:
  Q *point0, *point1;
  unsigned char size;
public:
  Ridge(Q *point0, Q *point1)
  : point0(point0), point1(point1), size(2)
  {}
  Ridge(Q *point0)
  : point0(point0), point1(NULL), size(1)
  {}
  Q* getPoint0() const {return point0;}
  Q* getPoint1() const {return point1;}
  unsigned char getSize() const {return size;}
};

class TerrainHoneycomb
{
private:
  std::vector<Q*> allQs;
  std::map<size_t, std::vector<Ridge>> cellRidges;
public:
  TerrainHoneycomb() = default;
  ~TerrainHoneycomb();
  void dumpQ(
    Point position, float elevation, size_t vorIndex,
    std::vector<size_t> nodes
  );
  void dumpNull();
  void dumpCellRidge(size_t cellID, Ridge ridge);
  Q* getQ(size_t idx);
  std::vector<Ridge> getCellRidges(size_t nodeID);
};

#endif