#ifndef TERRAINS_H
#define TERRAINS_H

#include <uchar.h>
#include <vector>

#include "point.hpp"

/**
 * @brief Represents a terrain primitive
 * 
 */
class T
{
private:
  Point loc;
  size_t cellID;
  float elevation;
public:
  /**
   * @brief Construct a new T object
   * 
   * @param loc Location
   * @param cellID ID of the hydrology primitive in which it is located
   */
  T(Point loc, size_t cellID)
  : loc(loc), cellID(cellID)
  {};
  /**
   * @brief Get the location of this primitive
   * 
   * @return const Point& The location
   */
  const Point& getLoc() const {return loc;}
  /**
   * @brief Get the ID of the hydrology primitive in which this primitive is located
   * 
   * @return size_t The ID of the hydrology primitive
   */
  size_t getCellID() const {return cellID;}
  /**
   * @brief Get the elevation of this primitive
   * 
   * @return float The elevation (in meters)
   */
  float getElevation() const {return elevation;}
  /**
   * @brief After the elevation of this primitive has been computed, use this method to set the new elevation
   * 
   * @param newElevation The elevation of this primitive
   */
  void setElevation(float newElevation) {elevation = newElevation;}
};

/**
 * @brief Analogous to the Terrain object in Python, but only a subset of the functionality. This class is basically just a wrapper for a vector or Primitives
 * 
 */
class Terrain
{
private:
  std::vector<T*> allTs;
public:
  ~Terrain();
  /**
   * @brief Create a terrain primitive and append it to the vector
   * 
   * This is intended for creating primitives from a binary data
   * stream
   * 
   * @param loc Location of the new primitive
   * @param cellID ID of the hydrology primitive in which it is located
   */
  void dumpT(Point loc, size_t cellID);
  /**
   * @brief The number of terrain primitives held by this object
   * 
   * @return size_t The number of primitives
   */
  size_t numTs() const;
  /**
   * @brief Gets the primitive identified by its index in the vector
   * 
   * @param idx The index of the primitive within the vector
   * @return T& A reference to the primitive
   */
  T& getT(size_t idx);
};

#endif