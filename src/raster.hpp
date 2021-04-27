#ifndef RASTER_H
#define RASTER_H

#include <iostream>
#include <math.h>

/*
2D array syntax is difficult in C++, especially when having to pass
arrays between functions. So this class basically emulates the
behavior of a 2D array.
*/
//template<typename T>
class Raster {
    private:
    size_t size;
    float resolution;
    float *data;

    public:
    Raster();
    ~Raster();
    void set(size_t col, size_t row, float datum);
    void set(float datum);
    void setResolution(float resolution);
    float get(float x, float y);
    void setSize(size_t size);
    size_t getRows();
    size_t getColumns();
};

#endif