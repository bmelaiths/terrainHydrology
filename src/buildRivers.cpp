#include <stdio.h>
#include <vector>
#include <endian.h>
#include <stdlib.h>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

#include "hydrologyParameters.hpp"
#include "hydrologyFunctions.hpp"

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

void openCVtest(HydrologyParameters& params) {
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

  cv::waitKey(50000);
}

int main() {
  //gather inputs

  #define INPUT stdin
  // #define FILEINPUT
  #ifdef FILEINPUT
  FILE *input = fopen("./binaryFile", "rb");

  if (input == NULL)
  {
    printf("Unable to open file\n");
    exit(1);
  }
  #endif

  HydrologyParameters params = readParamsFromStream(INPUT);

  // for(size_t candidate : params.candidates)
  // {
  //   printf
  //   (
  //     "Candidate %ld: (%f, %f)\n",
  //     candidate,
  //     params.hydrology.indexedNodes[candidate].loc.x,
  //     params.hydrology.indexedNodes[candidate].loc.y
  //   );
  // }

  /* Create a shoreline */
  // const int r = 1000;
  // cv::Mat src = cv::Mat::zeros( cv::Size( 4*r, 4*r ), CV_8U );
  // std::vector<cv::Point2f> vert(6);
  // vert[0] = cv::Point( 1500, 1340 );
  // vert[1] = cv::Point( 1000, 2000 );
  // vert[2] = cv::Point( 1500, 2860 );
  // vert[3] = cv::Point( 2500, 2860 );
  // vert[4] = cv::Point( 3000, 2000 );
  // vert[5] = cv::Point( 2500, 1340 );
  // for( int i = 0; i < 6; i++ )
  // {
  //   cv::line( src, vert[i],  vert[(i+1)%6], cv::Scalar( 255 ), 3 );
  // }
  // std::vector<std::vector<cv::Point> > contours;
  // cv::findContours( 
  //   src, contours, cv::RETR_TREE, cv::CHAIN_APPROX_SIMPLE
  // );

  // HydrologyParameters params;
  // params.Pa = 0.3;
  // params.Pc = 0.3;
  // params.maxTries = 15;
  // params.riverAngleDev = 1.0;
  // params.edgeLength = 10.0;
  // params.sigma = 0.75;
  // params.eta = 0.75;
  // params.zeta = 50.0;
  // params.slopeRate = 1.0;
  // params.resolution = 1.0;
  // params.riverSlope.setResolution(1.0);
  // params.riverSlope.setSize(3500);
  // params.riverSlope.set(1.0);
  // params.contour = contours[0];
  // const size_t contourIndex = 175;
  // params.candidates.push_back(
  //   params.hydrology.addMouthNode(
  //     Point(params.contour[contourIndex].x,params.contour[contourIndex].y),
  //     0.0, 1, contourIndex
  //   ).id
  // );

  // perform computatons

  // openCVtest(params);

  const uint8_t anotherNode = 0x2e, allDone = 0x21;

  while (
    params.candidates.size() > 0 // &&
    // params.hydrology.indexedNodes.size() < 100
  )
  {
    Primitive selectedCandidate = selectNode(params);
    alpha(selectedCandidate, params);
    // printf("\tNodes generated: %ld\r", params.hydrology.indexedNodes.size());
    fwrite(&anotherNode, sizeof(uint8_t), 1, stdout);
    fflush(stdout);
  }
  fwrite(&allDone, sizeof(uint8_t), 1, stdout);
  fflush(stdout);

  //export outputs

  writeBinary(params.hydrology, stdout);

  // printf("Number of nodes: %ld\n", params.hydrology.indexedNodes.size());


  //free resources
  #ifdef FILEINPUT
  fclose(INPUT);
  #endif

  return 0;
}