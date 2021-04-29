#ifndef HYDROPARAMS_H
#define HYDROPARAMS_H

#include <vector>
#include <random>

#include <opencv2/imgproc.hpp>

#include "raster.hpp"
#include "hydrology.hpp"

class HydrologyParameters
{
    public:
    float Pa, Pc;
    unsigned int maxTries;
    float riverAngleDev;
    float edgeLength;
    float sigma, eta, zeta;
    float slopeRate;
    float resolution;

    Raster<float> riverSlope;

    std::vector<cv::Point> contour;

    std::vector<Primitive*> candidates;
    Hydrology hydrology;

    std::default_random_engine generator;
    std::normal_distribution<float> distribution;
};

HydrologyParameters readParamsFromStream(FILE *stream);

float float_swap(float value);

#endif