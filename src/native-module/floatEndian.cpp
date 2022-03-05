#include "floatEndian.hpp"

#include <endian.h>
#include <stdint.h>

//this function written by user "synthetix" on cboard.cprogramming.com
float float_swap(float value) {
    union v {
        float f;
        uint32_t i;
    };

    union v val;

    val.f = value;
    val.i = be32toh(val.i);

    return val.f;
}

float float_tobe(float value) {
    union v {
        float f;
        uint32_t i;
    };

    union v val;

    val.f = value;
    val.i = htobe32(val.i);

    return val.f;
}