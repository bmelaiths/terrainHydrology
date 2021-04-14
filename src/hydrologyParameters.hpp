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

    Raster riverSlope;

    std::vector<Primitive> candidates;

    float resolution;

    std::vector<cv::Point> contour;

    public:
    HydrologyParameters(FILE *stream);
    
};

#endif