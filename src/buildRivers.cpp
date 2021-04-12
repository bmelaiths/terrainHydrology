#include <stdio.h>
#include <vector>

#include <opencv2/imgproc.hpp>
#include <opencv2/highgui.hpp>

/*
A list of data structures that will be used:
* A KDTree that will be used to query node locations
* A vector that will store actual node structs (this is where
  the node information will be, such as elevation, priority,
  etc). The tree will merely store an index in the vector
* Some data structure to represent the shore
*/

int main() {
    //gather inputs
    // char i, j;
    // fread(&i, 1, 1, stdin);
    // fread(&j, 1, 1, stdin);
    int64_t numPoints;
    int64_t xSize, ySize;
    fread(&numPoints, sizeof(int64_t), 1, stdin);
    fread(&xSize, sizeof(int64_t), 1, stdin);
    fread(&ySize, sizeof(int64_t), 1, stdin);

    #define X 0
    #define Y 1
    float inPoints[numPoints][2];
    for (int64_t i = 0; i < numPoints; i++)
    {
      // TODO adjust for endianness and size for differnt hardware
      // Points are recorded as y,x
      fread(inPoints[i]+Y, sizeof(float), 1, stdin);
      fread(inPoints[i]+X, sizeof(float), 1, stdin);
      // printf("Point %ld: (%f, %f)\n", i, inPoints[i][0], inPoints[i][1]);
    }


    //perform necessary computations
    // int result = i << j;
    std::vector<cv::Point> contour;
    for (int64_t i = 0; i < numPoints; i++)
    {
      // Points in a contour array are y,x
      contour.push_back(cv::Point2f(inPoints[i][X],inPoints[i][Y]));
    }

    cv::Mat raw_dist(cv::Size(xSize, ySize), CV_32F);
    for (int64_t y = 0; y < ySize; y++)
    {
      for (int64_t x = 0; x < xSize; x++)
      {
        raw_dist.at<float>(x,y) = (float) cv::pointPolygonTest(
          contour,
          cv::Point2f((float) x, (float) y),
          true
        );
      }
    }

    double minVal, maxVal;
    cv::Point maxDistPt; //inscribed circle center
    cv::minMaxLoc(raw_dist, &minVal, &maxVal, NULL, &maxDistPt);
    minVal = abs(minVal);
    maxVal = abs(maxVal);

    cv::Mat drawing = cv::Mat::zeros(cv::Size(xSize, ySize), CV_8UC3);
    for (int64_t y = 0; y < ySize; y++)
    {
      for (int64_t x = 0; x < xSize; x++)
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

    cv::waitKey(2000);

    //export outputs
    int result = 1;
    fwrite(&result, 1, 1, stdout);
    // printf("Number of points: %ld\n", numPoints);
    // printf("X Size: %ld\tY Size: %ld\n", xSize, ySize);
    // for (int64_t i = 0; i < numPoints; i++)
    // {
    //   printf("Point %ld: (%d, %d)\n", i, contour[i].x, contour[i].y);
    // }

    //free resources

    return 0;
}