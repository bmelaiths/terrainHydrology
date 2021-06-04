#ifndef TERRAIN_PRIMITIVES_H
#define TERRAIN_PRIMITIVES_H

#include <vector>

#include <opencv2/imgproc.hpp>

#include "hydrology.hpp"
#include "terrainHoneycomb.hpp"
#include "ts.hpp"

/**
 * @brief Reads data from the Python module
 * 
 */
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
  /**
   * @brief Construct a new Primitive Parameters object
   * 
   * @param stream A stream from which to receive the data
   */
  PrimitiveParameters(FILE *stream);
  ~PrimitiveParameters() = default;
};

#endif