#ifndef SHORE_H
#define SHORE_H

#include <vector>

#include <opencv2/imgproc.hpp>

class Shore {
    private:
        float resolution;
        size_t rows, cols;
        std::vector<cv::Point> contour;

        void toImageCoordinates(float x, float y, float *imageX, float *imageY);
    public:
        Shore();
        Shore(std::vector<cv::Point> contour, float resolution, size_t rasterXsize, size_t rasterYsize);

        double distanceToShore(float x, float y);
        cv::Point operator[](size_t idx);
};

#endif