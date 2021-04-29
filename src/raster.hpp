#ifndef RASTER_H
#define RASTER_H

#include <iostream>
#include <math.h>

/*
2D array syntax is difficult in C++, especially when having to pass
arrays between functions. So this class basically emulates the
behavior of a 2D array.
*/
template<typename T>
class Raster {
    private:
    size_t size;
    float resolution;
    T *data;

    public:
    Raster();
    ~Raster();
    void set(size_t col, size_t row, T datum);
    void set(T datum);
    void setResolution(float resolution);
    T get(float x, float y);
    void setSize(size_t size);
    size_t getRows();
    size_t getColumns();
};

template<typename T>
Raster<T>::Raster() : size(0), data(NULL) {}

template<typename T>
Raster<T>::~Raster() {
    if (data != NULL) {
        delete[] data;
    }
}

template<typename T>
void Raster<T>::set(size_t col, size_t row, T datum) {
    data[row * size + col] = datum;
}

template<typename T>
void Raster<T>::set(T datum) {
    for (size_t row = 0; row < size; row++)
    {
        for (size_t col = 0; col < size; col++)
        {
            set(col, row, datum);
        }
    }
}

template<typename T>
void Raster<T>::setResolution(float resolution) {
    this->resolution = resolution;
}

template<typename T>
T Raster<T>::get(float x, float y) {
    size_t col = (size_t) (x / resolution);
    size_t row = (size_t) (y / resolution);
    return data[row * size + col];
}

template<typename T>
void Raster<T>::setSize(size_t size) {
    this->size = size;
    try
    {
        data = new float[size * size];
    }
    catch (std::bad_alloc &ba)
    {
        std::cerr << "bad_alloc caught: " << ba.what();
    }
}

template<typename T>
size_t Raster<T>::getRows() {
    return size;
}

template<typename T>
size_t Raster<T>::getColumns() {
    return size;
}

#endif