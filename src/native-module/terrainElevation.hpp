#ifndef ELEVATION_H
#define ELEVATION_H

#include "hydrology.hpp"
#include "terrainHoneycomb.hpp"
#include "terrainPrimitives.hpp"
#include "point.hpp"
#include "shore.hpp"

/**
 * @file terrainElevation.hpp
 * 
 * The functions for computing elevations for terrain primitives
 */

/**
 * @brief Compute the elevation of a terrain primitive
 * 
 * @param t The primitive for which to compute
 * @param hydrology The Hydrology containing at least the relevant nodes
 * @param cells The TerrainHoneycomb containing at least the relevant ridges
 * @param ts The Terrain object (Not used)
 * @param contour The shoreline
 * @param resolution The number of meters per pixel in the raster that the contour is derived from
 * @param geosContext The GEOS context handle for the thread in which this function is called
 * @return float The elevation that the primitive should have
 */
float computePrimitiveElevation
(
  T& t, Hydrology& hydrology, TerrainHoneycomb& cells, Terrain& ts,
  Shore shore, float resolution,
  GEOSContextHandle_t geosContext
);

/**
 * @brief The distance to the line, and whether or nor the line is closest to an endpoint or the middle of the line
 * 
 */
struct endpointAndDistance
{
  float dist;
  bool isEndpoint;
};
/**
 * @brief Gets the distance between a point and a line
 * 
 * @param pLoc The point
 * @param l0Loc End 0 of the line
 * @param l1Loc End 1 of the line
 * @return endpointAndDistance An endpointAndDistance struct that contains not only the distance, but whether or not the closest point is one of the line's endpoints
 */
endpointAndDistance point_segment_distance
(const Point& pLoc, const Point& l0Loc, const Point& l1Loc);

/**
 * @brief Interpolates the elevation of the point on the ridge that a terrain primitive is closest to
 * 
 * This is a linear interpolation between q0 and q1
 * 
 * @param q0 End 0 of the ridge
 * @param q1 The other end of the ridge
 * @param t The terrain primitive
 * @param dist The distance between the terrain primitive and the ridge
 * @return float The elevation of the closest point (in meters)
 */
float lerpRidge(Q *q0, Q *q1, T& t, float dist);

/**
 * @brief The distance between two points
 * 
 * @param p0 
 * @param p1 
 * @return float 
 */
float distance(const Point& p0, const Point& p1);

#endif