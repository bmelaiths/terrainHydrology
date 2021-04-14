#include <stdio.h>
#include <vector>
#include <endian.h>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "point.hpp"
#include "kdtree.hpp"
#include "hydrology.hpp"
#include "hydrologyParameters.hpp"

/*
A list of data structures that will be used:
* A KDTree that will be used to query node locations
* A vector that will store actual node structs (this is where
  the node information will be, such as elevation, priority,
  etc). The tree will merely store an index in the vector
* Some data structure to represent the shore
* A graph
* A list of candidate nodes (this must be parallel)
*/

int main() {
  //gather inputs

  // FILE *input = fopen("./binaryFile", "rb");

  // if (input == NULL)
  // {
  //   printf("Unable to open file\n");
  //   exit(1);
  // }

  HydrologyParameters params(stdin);


  // perform computatons

  cv::Mat raw_dist(
    cv::Size(params.riverSlope.getColumns(), params.riverSlope.getRows()),
    CV_32F
  );

    #pragma omp parallel
    {
    #pragma omp for
    for (size_t y = 0; y < params.riverSlope.getRows(); y++)
    {
      for (size_t x = 0; x < params.riverSlope.getColumns(); x++)
      {
        raw_dist.at<float>(x,y) = (float) cv::pointPolygonTest(
          params.contour,
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

    cv::Mat drawing = cv::Mat::zeros(cv::Size(params.riverSlope.getColumns(), params.riverSlope.getRows()), CV_8UC3);
    for (size_t y = 0; y < params.riverSlope.getRows(); y++)
    {
      for (size_t x = 0; x < params.riverSlope.getColumns(); x++)
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


  //free resources

  // fclose(input);

    return 0;
}