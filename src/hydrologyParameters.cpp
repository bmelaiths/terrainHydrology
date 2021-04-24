#include "hydrologyParameters.hpp"

#include <stdio.h>
#include <endian.h>

//this function written by user "synthetix" on cboard.cprogramming.com
float float_swap(float value) {
    union v {
        float f;
        uint32_t i;
    };

    union v val;

    val.f = value;
    val.i = be32toh(val.i);

    return val.f;
}

HydrologyParameters readParamsFromStream(FILE *stream)
{
  HydrologyParameters params;

  fread(&params.Pa,               sizeof(float), 1, stream);
  fread(&params.Pc,               sizeof(float), 1, stream);
  fread(&params.edgeLength,       sizeof(float), 1, stream);
  fread(&params.sigma,            sizeof(float), 1, stream);
  fread(&params.eta,              sizeof(float), 1, stream);
  fread(&params.zeta,             sizeof(float), 1, stream);
  fread(&params.slopeRate,        sizeof(float), 1, stream);
  uint16_t maxTriesIn;
  fread(&maxTriesIn,      sizeof(uint16_t), 1, stream);
  fread(&params.riverAngleDev,    sizeof(float), 1, stream);
  params.Pa =            float_swap(params.Pa);
  params.Pc =            float_swap(params.Pc);
  params.edgeLength =    float_swap(params.edgeLength);
  params.sigma =         float_swap(params.sigma);
  params.eta =           float_swap(params.eta);
  params.zeta =          float_swap(params.zeta);
  params.slopeRate =     float_swap(params.slopeRate);
  params.maxTries =   (int) be16toh(maxTriesIn);
  params.riverAngleDev = float_swap(params.riverAngleDev);

  // printf("Pa: %f\n", Pa);
  // printf("Pc: %f\n", Pc);
  // printf("Edge Length: %f\n", edgeLength);
  // printf("Sigma: %f\n", sigma);
  // printf("Eta: %f\n", eta);
  // printf("Slope Rate: %f\n", slopeRate);
  // printf("Max Tries: %d\n", maxTries);
  // printf("River Angle Deviation: %f\n", riverAngleDev);

  uint32_t rasterXsize, rasterYsize;
  fread(&rasterXsize, sizeof(uint32_t), 1, stream);
  fread(&rasterYsize, sizeof(uint32_t), 1, stream);
  rasterXsize = be32toh(rasterXsize);
  rasterYsize = be32toh(rasterYsize);
  float riverSlopeIn;
  params.riverSlope.setSize(rasterXsize);
  for (uint32_t y = 0; y < rasterYsize; y++)
  {
    for (uint32_t x = 0; x < rasterXsize; x++)
    {
      fread(&riverSlopeIn, sizeof(float), 1, stream);
      params.riverSlope.set(x, y, float_swap(riverSlopeIn));
      // printf("Point (%d, %d): %f\n", x, y, riverSlopeIn);
    }
  }

  uint32_t numCandidates;
  fread(&numCandidates, sizeof(uint32_t), 1, stream);
  numCandidates = be32toh(numCandidates);
  // printf("Number of river mouths: %u\n", numCandidates);
  for (uint32_t i = 0; i < numCandidates; i++)
  {
    float x, y;
    uint32_t priority;
    uint64_t contourIndex;

    fread(&x, sizeof(float), 1, stream);
    fread(&y, sizeof(float), 1, stream);
    fread(&priority, sizeof(uint32_t), 1, stream);
    fread(&contourIndex, sizeof(uint64_t), 1, stream);

    x = float_swap(x);
    y = float_swap(y);
    priority = be32toh(priority);
    contourIndex = be64toh(contourIndex);

    params.candidates.push_back(
      params.hydrology.addMouthNode(
        Point(x,y), 0.0f, priority, contourIndex
      )
    );

    // printf("Node %d: (%f, %f), priority %d, contour index %ld\n",
    //     i, x, y, priority, contourIndex
    // );
  }

  uint64_t contourLength;
  fread(&params.resolution, sizeof(float), 1, stream);
  fread(&contourLength, sizeof(uint64_t), 1, stream);
  params.resolution = float_swap(params.resolution);
  contourLength = be64toh(contourLength);
  // printf("Spatial resolution: %f\n", resolution);
  // printf("Number of contour points: %ld\n", contourLength);
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
    params.contour.push_back(cv::Point2f(inPoints[i][X],inPoints[i][Y]));
  }

  params.riverSlope.setResolution(params.resolution);

  return params;
}