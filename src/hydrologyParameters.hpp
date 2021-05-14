#ifndef HYDROPARAMS_H
#define HYDROPARAMS_H

#include <vector>
#include <random>

#include <opencv2/imgproc.hpp>

#include "raster.hpp"
#include "hydrology.hpp"

/**
 * @brief A struct that holds all the necessary parameters to generate the river network
 * 
 */
class HydrologyParameters
{
    private:
    omp_lock_t candidateVectorLock;

    public:
    HydrologyParameters();
    /**
     * @brief Construct a new Hydrology Parameters object
     * 
     * @param lowerLeft The lower left corner of the expected area
     * @param upperRight The upper right corner of the expected area
     */
    HydrologyParameters(Point lowerLeft, Point upperRight);
    /**
     * @brief Construct a new Hydrology Parameters object with binary parameters coming from a stream
     * 
     * The stream should feed in all the parameters, including the
     * raster data.
     * 
     * @param stream A FILE stream, either stdin or a file containing the binary data
     */
    HydrologyParameters(FILE *stream);
    ~HydrologyParameters();
    HydrologyParameters(const HydrologyParameters& other);
    HydrologyParameters(HydrologyParameters&& other);

    HydrologyParameters& operator=(const HydrologyParameters& other);
    HydrologyParameters& operator=(HydrologyParameters&& other);

    float Pa, Pc;
    unsigned int maxTries;
    float riverAngleDev;
    float edgeLength;
    float sigma, eta, zeta;
    float slopeRate;
    float resolution;

    Raster<float> riverSlope;

    std::vector<cv::Point> contour;

    /**
     * @brief Acquires a lock on the vector of candidate nodes
     */
    void lockCandidateVector();
    std::vector<Primitive*> candidates;
    /**
     * @brief Releases the lock on the vector of candidate nodes
     */
    void unlockCandidateVector();
    Hydrology hydrology;

    std::default_random_engine generator;
    std::normal_distribution<float> distribution;
};

/**
 * @brief Converts a float from network order to system order
 */
float float_swap(float value);

#endif