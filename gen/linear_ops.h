#ifndef LINEAR_OPS_H_
#define LINEAR_OPS_H_
/**
 * \file linear_ops.h
 * \brief Linear encoding/decoding using matrix operations
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */
#include <stdint.h>
#include "gf2matrix.h"
/**
 * \brief
 * create and return a block-invertiable matrix of given size (bits)
 * and its inverse.
 * According to the algorithm in:
 * "Generating Large Non-Singular Matrices over an Arbitrary Field ..."
 * by Xiao, et al
 */
int make_block_invertible_matrix_pair(gf2matrix **m, gf2matrix **minv, int bits);
/**
 * \brief
 * create and return a block-invertible matrix of a given size
 * and initialize each block to Identity
 */
gf2matrix *make_block_identity_matrix(int bits);

/**
 * \brief
 * multiply a byte (8x1 vector) by a given 128x8
 * The result will be written out as a byte array starting
 * with the most signifcant bit
 */
int mul_byte_by_matrix_128x8(uint8_t result[16], const gf2matrix *m,
		uint8_t byte);
/**
 * \brief
 * multiply a 16-byte array (1x128 vector) by a given 128x128
 * The result (a 1x128 vector) will be written out as a byte array starting
 * with the most signifcant bit
 */
int mul_array_by_matrix_128x128(uint8_t result[16], const gf2matrix *m,
		const uint8_t value[16]);
/**
 * \brief
 * multiply a byte (8x1 vector) by a given 32x8
 * The result will be written out as a byte array starting
 * with the most signifcant bit
 */
int mul_byte_by_matrix_32x8(uint8_t *result, const gf2matrix *m,
		uint8_t byte);
/**
 * \brief
 * multiply a byte (8x1 vector) by a given nx8
 * where n is a multiple of 8.
 * The result will be written out as a byte array
 * of size n starting with the most signifcant bit
 */
int mul_byte_by_matrix(uint8_t *result, const gf2matrix *m,
		uint8_t byte);
/**
 * \brief
 * multiply a 4-byte array (1x32 vector) by a given 32x32
 * The result (a 1x32 vector) will be written out as a byte array starting
 * with the most signifcant bit
 */
int mul_array_by_matrix_32x32(uint8_t result[4], const gf2matrix *m,
		const uint8_t value[4]);
/**
 * \brief
 * Slice  src vertically into count submatrices (slices)
 * Destination slice matrices must have been fully allocated
 */
int slice_matrix_vertically(gf2matrix *slices[], int count,
		const gf2matrix *src);
/**
 * \brief
 * Allocate sliced matrices for a given  and slice count
 */
gf2matrix **new_sliced_vertical(int count, const gf2matrix *src);
/**
 * \brief de-allocate a given matrix array
 */
void free_matrices(gf2matrix **m, int count);
/**
 * \brief
 * wrapper for new_matrix
 */
gf2matrix *new_vector(int rows);
/**
 * \brief
 * wrapper for free_matrix
 */
void free_vector(gf2matrix *v);
/**
 * \brief
 * Copies bits from a vector (hi first) unto an integer
 * starting from the given offset row
 * Example [1 0 0 1 0 1 1] is the same as 0x0B or 11 decimal
 * if started with offset 3
 */
int vector2byte_offset(uint8_t *byte, const gf2matrix *v,
		int offset_row);
/**
 * \brief
 * Helper function to copy a vector v to an integer
 */
int vector2byte(uint8_t *byte, const gf2matrix *v);
/**
 * \brief
 * Helper function to copy the bits from a given integer to a vector v
 * Example 14 decimal (0x0E) will be copied as [1 1 1 0] to a vector
 */
int byte2vector_offset(gf2matrix *v, uint8_t byte, int offset_row);
/**
 * \brief
 * Helper function to copy the bits from a given integer to a vector v
 * Example 14 decimal (0x0E) will be copied as [1 1 1 0] to a vector
 */
int byte2vector(gf2matrix *v, uint8_t byte);

#endif /* LINEAR_OPS_H_ */
