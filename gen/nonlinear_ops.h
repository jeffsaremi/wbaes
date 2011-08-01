#ifndef NON_LINEAR_OPS_H_
#define NON_LINEAR_OPS_H_
/**
 * \file nonlinear_ops.h
 * \brief Non-linear (encoding/decoding using 4bit sboxes)
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */
#include <stdint.h>

/**
 * \brief
 * determines if a given bit vector satisifes the SAC
 * (Strict Avalanche Criterion)
 */
int is_bit_vector_sac(const uint8_t *f, int n);
/**
 * \brief
 * determines if a given bit vector satisifes the SAC
 * (Strict Avalanche Criterion)
 */
int is_sbox_sac(const uint8_t *f, int n);
/**
 * \brief
 * Extends a given Sbox F with a bit vector G which satisfies the SAC *
 * This method (and its notation) is based on:
 * "A Recursive Construction Method of Sboxes Satisfying Strict Avalanche Criterion"
 * by Kwangjo Kim, et al.
 */
uint8_t *extend_f_with_g(uint8_t *dest, const uint8_t *g,
		const uint8_t *f, int n, int k);
/**
 * brief
 * Extends a given sbox F using a bit vector of itself as G
 *
 */
uint8_t *extract_fk(uint8_t *dest, const uint8_t *f,
		int n, int k);
/**
 * brief
 * Extends a given sbox F by 1 bit
 *
 */
uint8_t *extend_f_by1(uint8_t *dest, const uint8_t *f,
		int n, int k);
/**
 * brief
 * inverts a given sbox F. If dest is NULL then an sbox is created
 *
 */
int invert_sbox(uint8_t *invsbox, const uint8_t *sbox, int n);
/**
 * brief
 * creates a random sbox and its invert for the given number of bits
 *
 * Bits must not be less than 3
 *
 */
int make_random_sbox(uint8_t *sbox, uint8_t *invsbox, int bits);

#endif /* NON_LINEAR_OPS_H_ */
