#ifndef TERRAIN_PRIMITIVES_H
#define TERRAIN_PRIMITIVES_H

#include <vector>

#include <opencv2/imgproc.hpp>

#include "hydrology.hpp"
#include "terrainHoneycomb.hpp"
#include "ts.hpp"

class PrimitiveParameters
{
public:
  float edgeLength;
  float resolution;
  std::vector<cv::Point> contour;
  Hydrology hydrology;
  TerrainHoneycomb cells;
  Terrain ts;
public:
  PrimitiveParameters(FILE *stream);
  ~PrimitiveParameters() = default;
};

#endif