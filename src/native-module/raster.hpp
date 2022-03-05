#ifndef RASTER_H
#define RASTER_H

#include <iostream>
#include <math.h>

/**
 * @brief Encapsulates raster data
 * 
 * @tparam T
 */
template<typename T>
class Raster {
    private:
    size_t rows, cols;
    float resolution;
    T *data;

    void toImageCoordinates(float x, float y, size_t *imageX, size_t *imageY);

    public:
    Raster();
    /**
     * @brief Construct a new Raster object
     * 
     * @param x Width of the raster data
     * @param y Height of the raster data
     * @param resolution Spatial resolution of the data: meters per unit
     */
    Raster(size_t x, size_t y, float resolution);
    ~Raster();
    Raster(const Raster& other);
    Raster(Raster&& other) noexcept;

    Raster& operator=(const Raster& other);
    Raster& operator=(Raster&& other) noexcept;

    /**
     * @brief Sets the value at a certain location
     * 
     * @param x X location in image coordinates
     * @param y Y location in image coordinates
     * @param datum The value to be set at that location
     */
    void set(size_t x, size_t y, T datum);
    /**
     * @brief Sets a value over the entire raster
     * 
     * @param datum The value
     */
    void set(T datum);

    /**
     * @brief Gets the value at a certain _spatial_ location
     * 
     * @param x The x location, _in meters_
     * @param y The y location, _in meters_
     * @return T The value at that location
     */
    T get(float x, float y);

    /**
     * @brief Get the number of vertical units
     * 
     * @return size_t 
     */
    size_t getRows();
    /**
     * @brief Get the number of horzontal units
     * 
     * @return size_t 
     */
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

    rows = std::move(other.rows);
    cols = std::move(other.cols);
    resolution = std::move(other.resolution);
    data = std::move(other.data);

    other.rows = 0;
    other.cols = 0;
    other.resolution = 0;
    other.data = NULL;

    return *this;
}

template<typename T>
void Raster<T>::toImageCoordinates(float x, float y, size_t *imageX, size_t *imageY) {
    x /= resolution;
    x += cols * 0.5;
    *imageX = x;

    y /= resolution;
    y = rows * 0.5 - y;
    *imageY = y;
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
    size_t col, row;

    toImageCoordinates(x, y, &col, &row);

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