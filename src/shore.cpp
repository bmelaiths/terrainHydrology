#include "shore.hpp"

void Shore::toImageCoordinates(float x, float y, float *imageX, float *imageY) {
    x /= resolution;
    x += cols * 0.5;
    *imageX = x;

    y /= resolution;
    y = rows * 0.5 - y;
    *imageY = y;
}

Shore::Shore() { }

Shore::Shore(std::vector<cv::Point> contour, float resolution, size_t rasterXsize, size_t rasterYsize):
    resolution(resolution), rows(rasterYsize), cols(rasterXsize), contour(contour)
{}

double Shore::distanceToShore(float x, float y) {
    float imageX, imageY;
    toImageCoordinates(x, y, &imageX, &imageY);

    return resolution * cv::pointPolygonTest(
        contour,
        cv::Point2f(imageX, imageY),
        true
    );
}

Point Shore::operator[](size_t idx) {
    return Point(contour[idx].x,contour[idx].y);
}