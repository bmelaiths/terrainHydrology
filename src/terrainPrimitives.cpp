#include "terrainPrimitives.hpp"

#include <geos_c.h>

#include "floatEndian.hpp"

PrimitiveParameters::PrimitiveParameters(FILE *stream, GEOSContextHandle_t geosContext)
{
  /*
    This method is very similar to HydrologyParameters::HydrologyParameters()

    Though verbose, this method is very simple. Basically,
    1. A variable is declared (unless it is declared in the class)
    2. fread() data into that variable
    3. Adjust the endian order, if necessary
  */

  float minX, maxX, minY, maxY;
  fread(&minX, sizeof(float), 1, stream);
  fread(&minY, sizeof(float), 1, stream);
  fread(&maxX, sizeof(float), 1, stream);
  fread(&maxY, sizeof(float), 1, stream);
  minX = float_swap(minX);
  minY = float_swap(minY);
  maxX = float_swap(maxX);
  maxY = float_swap(maxY);

  fread(&edgeLength, sizeof(float), 1, stream);
  edgeLength = float_swap(edgeLength);

  fread(&resolution, sizeof(float), 1, stream);
  resolution = float_swap(resolution);

  hydrology = Hydrology(Point(minX,minY), Point(maxX,maxY), edgeLength);

  /* Read the contour vector, and encode it in a structure
     that OpenCV will understand
   */
  uint64_t contourLength;
  fread(&contourLength, sizeof(uint64_t), 1, stream);
  contourLength = be64toh(contourLength);
  uint64_t inPoints[contourLength][2];
  fread(inPoints, sizeof(uint64_t), contourLength * 2, stream);
  #define X 0
  #define Y 1
  for (uint64_t i = 0; i < contourLength; i++)
  {
    inPoints[i][Y] = be64toh(inPoints[i][Y]);
    inPoints[i][X] = be64toh(inPoints[i][X]);
  }
  for (uint64_t i = 0; i < contourLength; i++)
  {
    // Points in a contour array are y,x
    contour.push_back(cv::Point2f(inPoints[i][Y],inPoints[i][X]));
  }

  /*
    Read in the hydrology nodes
  */
  uint64_t numNodes;
  fread(&numNodes, sizeof(uint64_t), 1, stream);
  numNodes = be64toh(numNodes);
  for (uint64_t i = 0; i < numNodes; i++)
  {
    uint32_t nodeID;
    float x, y, elevation;
    uint32_t parentID;
    uint32_t contourIndex;

    fread(&nodeID, sizeof(uint32_t), 1, stream);
    fread(&x, sizeof(float), 1, stream);
    fread(&y, sizeof(float), 1, stream);
    fread(&elevation, sizeof(float), 1, stream);
    fread(&parentID, sizeof(uint32_t), 1, stream);
    fread(&contourIndex, sizeof(uint32_t), 1, stream);
    nodeID = be32toh(nodeID);
    x = float_swap(x);
    y = float_swap(y);
    elevation = float_swap(elevation);
    parentID = be32toh(parentID);
    contourIndex = be32toh(contourIndex);

    /*
      Read in the rivers, using the GEOS library
    */

    uint8_t numRivers;
    fread(&numRivers, sizeof(uint8_t), 1, stream);

    std::vector<GEOSGeometry*> rivers;

    for (uint8_t river = 0; river < numRivers; river++)
    {
      uint16_t numPoints;
      fread(&numPoints, sizeof(uint16_t), 1, stream);
      numPoints = be16toh(numPoints);

      GEOSCoordSequence *string = GEOSCoordSeq_create_r(geosContext, numPoints, 3);

      for (uint16_t point = 0; point < numPoints; point++)
      {
        float x,y,z;

        fread(&x, sizeof(float), 1, stream);
        fread(&y, sizeof(float), 1, stream);
        fread(&z, sizeof(float), 1, stream);
        x = float_swap(x);
        y = float_swap(y);
        z = float_swap(z);

        GEOSCoordSeq_setXYZ_r(geosContext, string, point, x, y, z);
      }

      rivers.push_back(GEOSGeom_createLineString_r(geosContext, string));
    }

    float localWatershed, inheritedWatershed, flow;

    fread(&localWatershed, sizeof(float), 1, stream);
    fread(&inheritedWatershed, sizeof(float), 1, stream);
    fread(&flow, sizeof(float), 1, stream);
    localWatershed = float_swap(localWatershed);
    inheritedWatershed = float_swap(inheritedWatershed);
    flow = float_swap(flow);

    if (parentID == nodeID)
    {
      hydrology.dumpMouthNode(
        Point(x,y), elevation, 0, contourIndex, rivers,
        inheritedWatershed, localWatershed, flow
      );
    }
    else
    {
      hydrology.dumpRegularNode(
        Point(x,y), elevation, 0, parentID, rivers,
        inheritedWatershed, localWatershed, flow
      );
    }
  }

  /*
    Read the ridge primitives
  */
  uint64_t numQs;
  fread(&numQs, sizeof(uint64_t), 1, stream);
  numQs = be64toh(numQs);
  for (uint64_t i = 0; i < numQs; i++)
  {
    uint8_t hasQ;
    fread(&hasQ, sizeof(uint8_t), 1, stream);
    if (!hasQ)
    {
      cells.dumpNull();
      continue;
    }
    
    float x, y, elevation;
    uint64_t vorIndex;

    fread(&x, sizeof(float), 1, stream);
    fread(&y, sizeof(float), 1, stream);
    fread(&elevation, sizeof(float), 1, stream);
    fread(&vorIndex, sizeof(uint64_t), 1, stream);
    x = float_swap(x);
    y = float_swap(y);
    elevation = float_swap(elevation);
    vorIndex = be64toh(vorIndex);

    uint8_t numNeighbors;
    fread(&numNeighbors, sizeof(uint8_t), 1, stream);

    std::vector<uint64_t> neighbors;
    for (uint8_t neighbor = 0; neighbor < numNeighbors; neighbor++)
    {
      uint64_t neighborID;
      fread(&neighborID, sizeof(uint64_t), 1, stream);
      neighborID = be64toh(neighborID);
      neighbors.push_back(neighborID);
    }

    cells.dumpQ(Point(x,y), elevation, vorIndex, neighbors);
  }

  /*
    Read in the cells ridges dictionary
  */

  uint64_t cellsToProcess;
  fread(&cellsToProcess, sizeof(uint64_t), 1, stream);
  cellsToProcess = be64toh(cellsToProcess);
  for (uint64_t i = 0; i < cellsToProcess; i++)
  {
    uint64_t cellID;
    fread(&cellID, sizeof(uint64_t), 1, stream);
    cellID = be64toh(cellID);
    // printf("Node ID: %ld, ", cellID);

    uint8_t numRidges;
    fread(&numRidges, sizeof(uint8_t), 1, stream);
    for (uint8_t ri = 0; ri < numRidges; ri++)
    {
      uint8_t verticesInRidge;
      fread(&verticesInRidge, sizeof(uint8_t), 1, stream);
      
      if (verticesInRidge > 1)
      {
        uint64_t r0, r1;
        
        fread(&r0, sizeof(uint64_t), 1, stream);
        fread(&r1, sizeof(uint64_t), 1, stream);
        r0 = be64toh(r0);
        r1 = be64toh(r1);

        Ridge ridge(cells.getQ(r0), cells.getQ(r1));
        cells.dumpCellRidge(cellID, ridge);
        // printf("%ld, ", r1);
      }
      else
      {
        uint64_t r0;
        
        fread(&r0, sizeof(uint64_t), 1, stream);
        r0 = be64toh(r0);

        Ridge ridge(cells.getQ(r0));
        cells.dumpCellRidge(cellID, ridge);
      }
    }
  }

  /*
    Read in the terrain primitives
  */

  uint64_t tsToProcess;
  fread(&tsToProcess, sizeof(uint64_t), 1, stream);
  tsToProcess = be64toh(tsToProcess);
  for (uint64_t i = 0; i < tsToProcess; i++)
  {
    float x, y;
    uint32_t cellID;

    fread(&x, sizeof(float), 1, stream);
    fread(&y, sizeof(float), 1, stream);
    fread(&cellID, sizeof(uint32_t), 1, stream);
    x = float_swap(x);
    y = float_swap(y);
    cellID = be32toh(cellID);

    ts.dumpT(Point(x,y), cellID);
  }
}