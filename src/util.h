#ifndef _UTIL_H_
#define _UTIL_H_
/**
 * \file util.h
 * \brief Utility and helper functions
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */

#include <stdlib.h>
#include <stdint.h>

/**
 * \brief
 * returns a random number in the inclusive range: [min,max]
 */
int get_random(int min, int max);
/**
 * \brief
 * Utility to print bytes in hex
 */
void dump_hex(const char *label, const uint8_t *buf, int size);
/**
 * \brief
 * Utility to print bytes in hex
 */
void dump_hex2(const char *label, const uint8_t *buf, int size);
/**
 * \brief
 * Utility to print bytes in the state matrix
 */
void dump_state(const char *label, uint8_t state[4][4]);
/**
 * \brief
 * Utility to print 4 bytes in hex
 */
void dump_word(const uint8_t w[4], const char *label);
/**
 * \brief
 * Utility to print bytes in decimal
 */
void dump_decimal(const char *label, const uint8_t *buf, int size);
/**
 * \brief
 * Utility to print 4bit strips of 4 bytes each
 */
void dump_4bit_strip32(const char *label, uint8_t _4bitStrip[4][4][4]);
/**
 * \brief
 * Utility to print 4bit strips of 16 bytes each
 */
void dump_4bit_strip128(const char *label, uint8_t _4bitStrip[4][4][16]);

#endif /* _UTIL_H_ */
