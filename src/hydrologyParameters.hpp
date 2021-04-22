#ifndef HYDROPARAMS_H
#define HYDROPARAMS_H

#include <vector>

#include <opencv2/imgproc.hpp>

#include "raster.hpp"
#include "hydrology.hpp"

class HydrologyParameters
{
    public:
    float Pa, Pc;
    int maxTries;
    float riverAngleDev;
    float edgeLength;
    float sigma, eta, zeta;
    float slopeRate;
    float resolution;

    Raster riverSlope;

    std::vector<cv::Point> contour;

    std::vector<Primitive> candidates;
    Hydrology hydrology;
};

HydrologyParameters readParamsFromStream(FILE *stream);

#endif