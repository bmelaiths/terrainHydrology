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

float Raster::get(size_t col, size_t row) {
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