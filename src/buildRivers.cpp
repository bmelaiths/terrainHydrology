#include <stdio.h>
#include <vector>
#include <endian.h>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "point.hpp"
#include "kdtree.hpp"
#include "hydrology.hpp"

/*
A list of data structures that will be used:
* A KDTree that will be used to query node locations
* A vector that will store actual node structs (this is where
  the node information will be, such as elevation, priority,
  etc). The tree will merely store an index in the vector
* Some data structure to represent the shore
* A graph
*/

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

int main() {
    //gather inputs

    float Pa, Pc;
    float edgeLength;
    float sigma, eta;
    float slopeRate;
    uint16_t maxTries;
    float riverAngleDev;

    fread(&Pa,               sizeof(float), 1, stdin);
    fread(&Pc,               sizeof(float), 1, stdin);
    fread(&edgeLength,       sizeof(float), 1, stdin);
    fread(&sigma,            sizeof(float), 1, stdin);
    fread(&eta,              sizeof(float), 1, stdin);
    fread(&slopeRate,        sizeof(float), 1, stdin);
    fread(&maxTries,      sizeof(uint16_t), 1, stdin);
    fread(&riverAngleDev,    sizeof(float), 1, stdin);
    Pa =            float_swap(Pa);
    Pc =            float_swap(Pc);
    edgeLength =    float_swap(edgeLength);
    sigma =         float_swap(sigma);
    eta =           float_swap(eta);
    slopeRate =     float_swap(slopeRate);
    maxTries =         be16toh(maxTries);
    riverAngleDev = float_swap(riverAngleDev);

    // printf("Pa: %f\n", Pa);
    // printf("Pc: %f\n", Pc);
    // printf("Edge Length: %f\n", edgeLength);
    // printf("Sigma: %f\n", sigma);
    // printf("Eta: %f\n", eta);
    // printf("Slope Rate: %f\n", slopeRate);
    // printf("Max Tries: %d\n", maxTries);
    // printf("River Angle Deviation: %f\n", riverAngleDev);

    uint32_t rasterXsize, rasterYsize;
    fread(&rasterXsize, sizeof(uint32_t), 1, stdin);
    fread(&rasterYsize, sizeof(uint32_t), 1, stdin);
    rasterXsize = be32toh(rasterXsize);
    rasterYsize = be32toh(rasterYsize);
    float riverSlope[rasterYsize][rasterXsize];
    fread(riverSlope, sizeof(float), rasterXsize * rasterYsize, stdin);

    // printf("Raster X Size: %d\n", rasterXsize);
    // printf("Raster Y Size: %d\n", rasterYsize);
    // uint32_t maxX = 0, maxY = 0;
    for (uint32_t y = 0; y < rasterYsize; y++)
    {
        for (uint32_t x = 0; x < rasterXsize; x++)
        {
            riverSlope[y][x] = float_swap(riverSlope[y][x]);
            // if (riverSlope[y][x] > riverSlope[maxY][maxX])
            // {
            //     maxX = x;
            //     maxY = y;
            // }
            // // printf("Point (%d, %d): %f\n", x, y, riverSlope[y][x]);
        }
    }
    // // printf("Max Value: %f\n", riverSlope[maxY][maxX]);

    uint32_t numCandidates;
    fread(&numCandidates, sizeof(uint32_t), 1, stdin);
    numCandidates = be32toh(numCandidates);
    // printf("Number of river mouths: %u\n", numCandidates);
    Primitive candidates[numCandidates];
    for (uint32_t i = 0; i < numCandidates; i++)
    {
        float x, y;
        uint32_t priority;
        uint64_t contourIndex;

        fread(&x, sizeof(float), 1, stdin);
        fread(&y, sizeof(float), 1, stdin);
        fread(&priority, sizeof(uint32_t), 1, stdin);
        fread(&contourIndex, sizeof(uint64_t), 1, stdin);

        x = float_swap(x);
        y = float_swap(y);
        priority = be32toh(priority);
        contourIndex = be64toh(contourIndex);

        candidates[i] = Primitive(Point(x,y), priority, contourIndex);

        // printf("Node %d: (%f, %f), priority %d, contour index %ld\n",
        //     i, x, y, priority, contourIndex
        // );
    }

    float resolution;
    uint64_t contourLength;
    fread(&resolution, sizeof(float), 1, stdin);
    fread(&contourLength, sizeof(uint64_t), 1, stdin);
    resolution = float_swap(resolution);
    contourLength = be64toh(contourLength);
    // printf("Spatial resolution: %f\n", resolution);
    // printf("Number of contour points: %ld\n", contourLength);
    uint64_t inPoints[contourLength][2];
    fread(inPoints, sizeof(uint64_t), contourLength * 2, stdin);
    #define X 0
    #define Y 1
    for (uint64_t i = 0; i < contourLength; i++)
    {
      inPoints[i][Y] = be64toh(inPoints[i][Y]);
      inPoints[i][X] = be64toh(inPoints[i][X]);
    }


    //perform necessary computations
    // int result = i << j;
    std::vector<cv::Point> contour;
    for (uint64_t i = 0; i < contourLength; i++)
    {
      // Points in a contour array are y,x
      contour.push_back(cv::Point2f(inPoints[i][X],inPoints[i][Y]));
    }

    cv::Mat raw_dist(cv::Size(rasterXsize, rasterYsize), CV_32F);

    #pragma omp parallel
    {
    #pragma omp for
    for (int64_t y = 0; y < rasterYsize; y++)
    {
      for (int64_t x = 0; x < rasterXsize; x++)
      {
        raw_dist.at<float>(x,y) = (float) cv::pointPolygonTest(
          contour,
          cv::Point2f((float) x, (float) y),
          true
        );
      }
    }
    }

    double minVal, maxVal;
    cv::Point maxDistPt; //inscribed circle center
    cv::minMaxLoc(raw_dist, &minVal, &maxVal, NULL, &maxDistPt);
    minVal = abs(minVal);
    maxVal = abs(maxVal);

    cv::Mat drawing = cv::Mat::zeros(cv::Size(rasterXsize, rasterYsize), CV_8UC3);
    for (int64_t y = 0; y < rasterYsize; y++)
    {
      for (int64_t x = 0; x < rasterXsize; x++)
      {
        if( raw_dist.at<float>(x,y) < 0 )
        {
          drawing.at<cv::Vec3b>(x,y)[0] = (uchar)(255 - abs(raw_dist.at<float>(x,y)) * 255 / minVal);
        }
        else if( raw_dist.at<float>(x,y) > 0 )
        {
          drawing.at<cv::Vec3b>(x,y)[2] = (uchar)(255 - raw_dist.at<float>(x,y) * 255 / maxVal);
        }
        else
        {
          drawing.at<cv::Vec3b>(x,y)[0] = 255;
          drawing.at<cv::Vec3b>(x,y)[1] = 255;
          drawing.at<cv::Vec3b>(x,y)[2] = 255;
        }
      }
    }

    cv::imshow("See if this works", drawing);

    cv::waitKey(10000);

    //export outputs
    uint8_t result = 1;
    fwrite(&result, sizeof(uint8_t), 1, stdout);
    // printf("Number of points: %ld\n", numPoints);
    // printf("X Size: %ld\tY Size: %ld\n", xSize, ySize);
    // for (int64_t i = 0; i < numPoints; i++)
    // {
      // printf("Point %ld: (%d, %d)\n", i, contour[i].x, contour[i].y);
    // }

    //free resources

    return 0;
}