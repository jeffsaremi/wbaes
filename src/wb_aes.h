#ifndef _WB_AES_H
#define _WB_AES_H
/**
 * \file wb_aes.h
 * \brief WhiteBox AES runtime (encryption/decryption) functions
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */

#include <stdint.h>
#include "aes.h"
#include "gf2matrix.h"

/** Type I Lookup tables */
typedef uint8_t (typeI_t)[4][4][256][16];
/** Type IA Lookup tables (same as typeI_t) */
typedef typeI_t typeIA_t;
/** Type IB Lookup tables (same as typeI_t) */
typedef typeI_t typeIB_t;
/** Type II or Type III Lookup tables (for one round) */
typedef uint8_t (typeII_round_t)[4][4][256][4];
/** Type II or Type III Lookup tables (for one round) */
typedef typeII_round_t typeIII_round_t;
/** Type II or Type III Lookup tables for all rounds */
typedef typeII_round_t (typeII_t)[NR-1];
/** Type II or Type III Lookup tables for all rounds */
typedef typeII_t typeIII_t;
/** Type IV Lookup tables following Type I (128 bit) */
typedef uint8_t (typeIV128_t)[4][4][15][2][128];
/** Type IV Lookup tables following Type IA (128 bit) */
typedef typeIV128_t typeIV_IA_t;
/** Type IV Lookup tables following Type IB (128 bit) */
typedef typeIV128_t typeIV_IB_t;
/** Type IV Lookup tables following Type II or III (32 bit) for one round */
typedef uint8_t (typeIV32_round_t)[4][4][3][2][128];
/** Type IV Lookup tables following Type II or III (32 bit) for all rounds */
typedef typeIV32_round_t (typeIV32_t)[NR-1];
/** Type IV Lookup tables following Type II (32 bit) for all rounds */
typedef typeIV32_round_t typeIV_II_round_t;
/** Type IV Lookup tables following Type III (32 bit) for all rounds */
typedef typeIV32_round_t typeIV_III_round_t;
/** 4-bit strips for 32bit operations */
typedef uint8_t _4bit_strip32_t[4][4][4];
/** 4-bit strips for 128bit operations */
typedef uint8_t _4bit_strip128_t[4][4][16];
/** 2 4x4sboxes used in 4bit encoding/decoding */
typedef uint8_t (sboxes_8bit_t)[4][4][2][16];

/**
 * \brief
 * substitute hi and low 4-bits in x using the given SBOXes
 */
uint8_t sub_bytes_hi_lo(uint8_t x, uint8_t _sbox[][16]);

/**
 * \brief
 * Perform XOR of a and b by looking up the result in the
 * given lookup table.
 * The arguments a and b are supposed to be 4bits each.
 */
uint8_t do_xor(uint8_t a, uint8_t b, uint8_t lookup_8x4[128]);
/**
 * \brief
 * Perform simulataneous XOR of two sets of 4bits (lo and hi bits).
 * Combine the result into an 8bit and return.
 * Two separate lookup tables are needed.
 */
uint8_t do_xor_hi_lo(uint8_t a, uint8_t b, uint8_t dual_lookup_8x4[][128]);
/**
 * \brief
 * Perform Type IV on the 128bit strips and save the result in the state.
 * 4bit strips are 32 4-bit values for each row/col of the previous step.
 * These 4bits will be XOR'd by using the provided TypeIV tables.
 */
void do_typeIV128(uint8_t state[4][4], _4bit_strip128_t strips,
		typeIV128_t typeIV_Is);
/**
 * \brief
 * Perform Type IV on the 32bit strips and save the result in the state.
 * 4bit strips are 8 4-bit values for each row/col of the previous step.
 * These 4bits will be XOR'd by using the provided TypeIV tables.
 */
void do_typeIV32(uint8_t state[4][4], _4bit_strip32_t strips,
		typeIV32_round_t round_typeIV32);
/**
 * \brief
 * Perform a type I operation on the given state using the supplied
 * Type I lookup tables.
 * The results are stored as 16 bytes or 32 4-bit values for each state
 * value
 */
void do_typeI(_4bit_strip128_t strips, uint8_t state[4][4],
		typeI_t typeIs);
/**
 * \brief
 * Perform a type II operation on the given state using the supplied
 * Type II lookup tables of the given round.
 * The results are stored as 4 bytes or 8 4-bit values for each state
 * value
 */
void do_typeII(_4bit_strip32_t strips, uint8_t state[4][4],
		typeII_round_t round_typeIIs);

/**
 * \brief
 * A convenience wrapper for do_typeIV32.
 */
void do_typeIV_II(uint8_t state[4][4], _4bit_strip32_t strips,
		typeIV_II_round_t round_typeIV_IIs);
/**
 * \brief
 * A convenience wrapper for do_typeIV32.
 */
void do_typeIV_III(uint8_t state[4][4], _4bit_strip32_t strips,
		typeIV_III_round_t round_typeIV_IIIs);
/**
 * \brief
 * Perform a type III operation (same as Type II) on the given state
 * using the supplied Type III lookup tables of the given round.
 * The results are stored as 16 bytes or 32 4-bit values for each state
 * value
 */
void do_typeIII(_4bit_strip32_t strips, uint8_t state[4][4],
		typeIII_round_t round_typeIIIs);
/**
 * \brief
 * A convenience wrapper for do_typeI.
 */
void do_typeIA(_4bit_strip128_t strips, uint8_t state[4][4],
		typeIA_t typeIAs);
/**
 * \brief
 * A convenience wrapper for do_typeI.
 */
void do_typeIB(_4bit_strip128_t strips, uint8_t state[4][4],
		typeIB_t typeIBs);
/**
 * \brief
 * A convenience wrapper for do_typeIV128.
 */
void do_typeIV_IA(uint8_t state[4][4], _4bit_strip128_t strips,
		typeIV_IA_t typeIV_IAs);
/**
 * \brief
 * A convenience wrapper for do_typeIV128.
 */
void do_typeIV_IB(uint8_t state[4][4], _4bit_strip128_t strips,
		typeIV_IB_t typeIV_IBs);

/**
 * \brief
 * Perform input encoding using the given linear and non-linear tables
 */
void do_input(uint8_t state[4][4], const uint8_t input[KEY_SIZE],
		gf2matrix *linear_encoding, sboxes_8bit_t input_sboxes);
/**
 * \brief
 * Perform output decoding using the given linear and non-linear tables
 */
void do_output(uint8_t output[KEY_SIZE], uint8_t state[4][4],
		gf2matrix *linear_decoding, sboxes_8bit_t output_sboxes);

#endif /*  _WB_AES_H */
