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
  /**
   * @brief Construct a new Q object using data from the Python stream
   * 
   * @param position This primitive's location
   * @param elevation This primitive's elevation
   * @param vorIndex This index of this Q's vertex. (Not used)
   * @param nodes IDs of the hydrology primitives that this Q borders
   */
  Q
  (
    Point position, float elevation,
    size_t vorIndex, std::vector<size_t> nodes
  ):
  position(position), elevation(elevation), vorIndex(vorIndex),
  nodes(nodes)
  {}
  ~Q() = default;
  /**
   * @brief Get the location of this Q
   * 
   * @return const Point The location
   */
  const Point getPosition() const {return position;}
  /**
   * @brief Get the elevation of this Q
   * 
   * @return float The elevation (in meters)
   */
  float getElevation() const {return elevation;}
  /**
   * @brief Get the index of the voronoi vertex that this Q is based on
   * 
   * @return size_t The index
   */
  size_t getVorIndex() const {return vorIndex;}
  /**
   * @brief Get the IDs of the hydrology nodes that this Q borders
   * 
   * @return std::vector<size_t> 
   */
  std::vector<size_t> getNodes() const {return nodes;}
};

/**
 * @brief Represents a ridge, comprised of one or two Q primitives
 * 
 */
class Ridge
{
private:
  Q *point0, *point1;
  unsigned char size;
public:
  /**
   * @brief Construct a new Ridge object from two Q primitives
   * 
   * @param point0 
   * @param point1 
   */
  Ridge(Q *point0, Q *point1)
  : point0(point0), point1(point1), size(2)
  {}
  /**
   * @brief Construct a new Ridge object from just one Q primitive
   * 
   * @param point0 
   */
  Ridge(Q *point0)
  : point0(point0), point1(NULL), size(1)
  {}
  /**
   * @brief Get end 0 of the ridge
   * 
   * @return Q* 
   */
  Q* getPoint0() const {return point0;}
  /**
   * @brief Get the other end of the ridge, if it exists
   * 
   * @return Q* 
   */
  Q* getPoint1() const {return point1;}
  /**
   * @brief Get the number of primitives in this ridge
   * 
   * @return unsigned char Either 1 or 2
   */
  unsigned char getSize() const {return size;}
};

/**
 * @brief This class associates Q primitives with cells and ridges
 * 
 * This class is a subset of the functionality of its analogous
 * Python class.
 * 
 */
class TerrainHoneycomb
{
private:
  std::vector<Q*> allQs;
  std::map<size_t, std::vector<Ridge>> cellRidges;
public:
  TerrainHoneycomb() = default;
  ~TerrainHoneycomb();

  /**
   * @brief Creates a Q primitive with the specified properties and appends it ot the vector
   * 
   * This is intended for creating primitives from a binary data stream
   * 
   * @param position The position of the new Q
   * @param elevation The elevation of the new Q
   * @param vorIndex The voronoi index of the new Q's corresponding vertex
   * @param nodes The hydrology primitives that this Q borders
   */
  void dumpQ(
    Point position, float elevation, size_t vorIndex,
    std::vector<size_t> nodes
  );
  /**
   * @brief Appends a null pointer to the internal vector
   * 
   * This is important because every Q must be at a specific position in
   * the vector. The binary data from the Python module that encodes the
   * cellRidges dictionary makes reference to specific indices in that
   * vector/list.
   */
  void dumpNull();

  /**
   * @brief Associates a Ridge with the ID of a hydrology primitive
   * 
   * @param cellID The ID of the hydrology primitive that this ridge encloses
   * @param ridge The ridge
   */
  void dumpCellRidge(size_t cellID, Ridge ridge);

  /**
   * @brief Gets a pointer to the Q primitive at an index within the vector
   * 
   * @param idx The index within the vector
   * @return Q* This pointer may be NULL if that is what is found at the index
   */
  Q* getQ(size_t idx);
  /**
   * @brief Gets the ridges that enclose a hydrology cell
   * 
   * @param nodeID The ID of the cell
   * @return std::vector<Ridge> The ridges that enclose it
   */
  std::vector<Ridge> getCellRidges(size_t nodeID);
};

#endif