#ifndef TEST_UTIL_H_
#define TEST_UTIL_H_

#include "wb_aes_gen.h"
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

/**
 * \brief
 * copies one column from the state table to the given array
 */
void copy_state_col(uint8_t col_bytes[4], uint8_t state[4][4], int col);
/**
 * \brief
 * copies a given vector over a column in the state table
 */
void copy_col2state(uint8_t state[4][4], uint8_t col_bytes[4], int col);
/**
 * \brief
 * assigns random numbers to state table members
 */
void randomize_state(uint8_t state[4][4]);
/**
 * \brief
 * performs an inverse of MixColums to a given state (used in decipher)
 */
void inv_mix_columns(uint8_t state[4][4]);
/**
 * \brief
 * creates an Identity T-Box
 */
void make_identity_tbox(tbox_t tbox);
/**
 * \brief
 * creates a 4-bit SBox which basically does not substitution
 */
void make_identity_4bit_sbox(uint8_t sbox[16]);
/**
 * \brief
 * creates a set of Identity SBoxes for 8-bit operations
 */
void make_identity_sboxes8(sboxes_8bit_t sboxes);
/**
 * \brief
 * creates a set of Identity SBoxes for 32-bit operations
 */
void make_identity_sboxes32(sboxes_32bit_t sboxes);
/**
 * \brief
 * creates a set of Identity SBoxes for 128-bit operations
 */
void make_identity_sboxes128(sboxes_128bit_t sboxes);
/**
 * \brief
 * creates a set of Identity T-Box Mixing Bijections
 */
void make_identity_tbox_mixing_bijections(tbox_mixing_bijections_t tbox_mbs);
/**
 * \brief
 * copies random values to given key of pre-defined size
 */
void randomize_key(uint8_t *key);



#endif /* TEST_UTIL_H_ */
