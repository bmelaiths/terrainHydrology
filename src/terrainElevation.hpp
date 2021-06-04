#ifndef ELEVATION_H
#define ELEVATION_H

#include "hydrology.hpp"
#include "terrainHoneycomb.hpp"
#include "terrainPrimitives.hpp"
#include "point.hpp"

float computePrimitiveElevation
(
  T& t, Hydrology& hydrology, TerrainHoneycomb& cells, Terrain& ts,
  std::vector<cv::Point>& contour, float resolution
);

struct endpointAndDistance
{
  float dist;
  bool isEndpoint;
};
endpointAndDistance point_segment_distance
(const Point& tLoc, const Point& q0Loc, const Point& q1Loc);

float lerpRidge(Q *q0, Q *q1, T& t, float dist);

float distance(const Point& p0, const Point& p1);

#endif