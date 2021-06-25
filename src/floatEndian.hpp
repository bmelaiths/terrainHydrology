#ifndef FLOAT_ENDIAN_H
#define FLOAT_ENDIAN_H

/**
 * @brief Converts a float from network order to system order
 */
float float_swap(float value);

/**
 * @brief Converts a float from system order to network order
 */
float float_tobe(float value);

#endif