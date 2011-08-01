#ifndef TEST_UTIL_H_
#define TEST_UTIL_H_

#include <stdint.h>

/**
 * \brief
 * compares two given state vectors. Return 0 if equal
 */
int comp_states(uint8_t state1[4][4], uint8_t state2[4][4]);
/**
 * \brief
 * compares two given 4bit strips of 4 bytes each. Return 0 if equal
 */
int comp_4bit_strips32(uint8_t _4bitStrip1[4][4][4],
		uint8_t _4bitStrip2[4][4][4]);
/**
 * \brief
 * compares two given 4bit strips of 16 bytes each. Return 0 if equal
 */
int comp_4bit_strips128(uint8_t _4bitStrip1[4][4][16],
		uint8_t _4bitStrip2[4][4][16]);


#endif /* TEST_UTIL_H_ */
