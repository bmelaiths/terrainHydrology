#include "raster.hpp"

Raster::Raster() : size(0), data(NULL) {}

Raster::~Raster() {
    if (data != NULL) {
        delete[] data;
    }
}

void Raster::set(size_t col, size_t row, float datum) {
    data[row * size + col] = datum;
}

void Raster::set(float datum) {
    for (size_t row = 0; row < size; row++)
    {
        for (size_t col = 0; col < size; col++)
        {
            set(col, row, datum);
        }
    }
}

void Raster::setResolution(float resolution) {
    this->resolution = resolution;
}

float Raster::get(float x, float y) {
    size_t row = (size_t) (x / resolution);
    size_t col = (size_t) (y / resolution);
    return data[row * size + col];
}

void Raster::setSize(size_t size) {
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

size_t Raster::getRows() {
    return size;
}

size_t Raster::getColumns() {
    return size;
}