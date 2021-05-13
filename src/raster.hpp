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
    size_t rows, cols;
    float resolution;
    T *data;

    public:
    Raster();
    Raster(size_t x, size_t y, float resolution);
    ~Raster();
    Raster(const Raster& other);
    Raster(Raster&& other) noexcept;

    Raster& operator=(const Raster& other);
    Raster& operator=(Raster&& other) noexcept;

    void set(size_t col, size_t row, T datum);
    void set(T datum);

    T get(float x, float y);

    size_t getRows();
    size_t getColumns();
};

template<typename T>
Raster<T>::Raster() : rows(0), cols(0), data(NULL) {}

template<typename T>
Raster<T>::Raster(size_t rows, size_t cols, float resolution)
: rows(rows), cols(cols), resolution(resolution)
{
    data = new float[rows * cols];
}

template<typename T>
Raster<T>::~Raster() {
    if (data != NULL) {
        delete[] data;
    }
}

template<typename T>
Raster<T>::Raster(const Raster& other)
: rows(other.rows), cols(other.cols), resolution(other.resolution)
{
    data = new T[rows * cols];
    memcpy(data, other.data, sizeof(T) * rows * cols);
}

template<typename T>
Raster<T>::Raster(Raster&& other) noexcept
: rows(std::move(other.rows)), cols(std::move(other.cols)),
  resolution(std::move(other.resolution)), data(std::move(other.data))
{
    other.rows = 0;
    other.cols = 0;
    other.data = NULL;
}

template<typename T>
Raster<T>& Raster<T>::operator=(const Raster<T>& other)
{
    if (this == &other)
    {
        return *this;
    }

    rows = other.rows;
    cols = other.cols;
    resolution = other.resolution;

    data = new T[rows * cols];
    memcpy(data, other.data, sizeof(T) * rows * cols);

    return *this;
}

template<typename T>
Raster<T>& Raster<T>::operator=(Raster<T>&& other) noexcept
{
    if (this == &other)
    {
        return *this;
    }

    rows = std::move(rows);
    cols = std::move(other.cols);
    resolution = std::move(other.resolution);
    data = std::move(other.data);

    other.rows = 0;
    other.cols = 0;
    other.data = NULL;

    return *this;
}

template<typename T>
void Raster<T>::set(size_t col, size_t row, T datum) {
    data[row * cols + col] = datum;
}

template<typename T>
void Raster<T>::set(T datum) {
    for (size_t row = 0; row < rows; row++)
    {
        for (size_t col = 0; col < cols; col++)
        {
            set(col, row, datum);
        }
    }
}

template<typename T>
T Raster<T>::get(float x, float y) {
    size_t col = (size_t) (x / resolution);
    size_t row = (size_t) (y / resolution);
    return data[row * cols + col];
}

template<typename T>
size_t Raster<T>::getRows() {
    return rows;
}

template<typename T>
size_t Raster<T>::getColumns() {
    return cols;
}

#endif