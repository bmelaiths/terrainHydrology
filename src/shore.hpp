#ifndef SHORE_H
#define SHORE_H

#include <vector>

#include <opencv2/imgproc.hpp>

#include "point.hpp"

class Shore {
    private:
        float resolution;
        size_t rows, cols;
        std::vector<cv::Point> contour;

        void toImageCoordinates(float x, float y, float *imageX, float *imageY);
    public:
        Shore();
        /**
         * @brief Construct a new Shore object
         * 
         * @param contour The openCV contour that represents the shoreline
         * @param resolution The resolution of the original shore input image in meters per pixel
         * @param rasterXsize The width of the original shore image in pixels
         * @param rasterYsize The height of the original shore image in pixels
         */
        Shore(std::vector<cv::Point> contour, float resolution, size_t rasterXsize, size_t rasterYsize);

        /**
         * @brief Gets the distance between a point and the shore
         * 
         * @param x The location to test (meters)
         * @param y The location to test (meters)
         * @return double The distance between the test location and the shore (meters)
         */
        double distanceToShore(float x, float y);
        /**
         * @brief Gets a point on the shore by index
         * 
         * The shore can be thought of as an array of points. These points demarcate the boundary between land and sea, and can be indexed
         * 
         * @param idx The index
         * @return Point The index-th coordinate of the shoreline
         */
        Point operator[](size_t idx);
};

#endif