#ifndef _MATRIX_H_
#define _MATRIX_H_

/**
 * \file gf2matrix.h
 * \brief GF(2) Matrix operations
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */
#include <stdint.h>
#include "m4ri_gf2matrix_typedef.h"

/**
 * \brief allocate a new matrix
 */
gf2matrix *new_matrix(int rows, int cols);
/**
 * \brief de-allocate a given matrix
 */
void free_matrix(gf2matrix *m);
/**
 * \brief
 * get row count of a matrix
 */
int get_rows(const gf2matrix *m);
/**
 * \brief
 * get column count of a matrix
 */
int get_cols(const gf2matrix *m);
/**
 * \brief
 * create a new matrix as an identity matrix
 */
gf2matrix *make_identity_matrix(int bits);
/**
 * \brief
 * convert an existing matrix into an identity matrix
 */
void convert_to_identity(gf2matrix *m);
/**
 * \brief returns 1 if an inverse for the given matrix can be found
 */
int is_invertible(const gf2matrix *m);
/**
 * \brief
 * invert a given matrix; use the destination if given or create a new matrix
 */
gf2matrix *invert_matrix(gf2matrix *dest, const gf2matrix *src);
/**
 * \brief
 * create a copy of a given matrix
 */
gf2matrix *dup_matrix(const gf2matrix *m);
/**
 * \brief
 * add A and B to give C; If C is NULL then create and return a new one
 * C = A + B
 */
gf2matrix *add_matrices(gf2matrix *C, const gf2matrix *A, const gf2matrix *B);
/**
 * \brief
 * multiply A and B to give C; If C is NULL then create and return a new one
 * C = A B
 */
gf2matrix *mul_matrices(gf2matrix *C, const gf2matrix *A, const gf2matrix *B);
/**
 * \brief
 * compare matrices; if equal return 0
 */
int comp_matrices(const gf2matrix *a, const gf2matrix *b);
/**
 * \brief
 * clear all bits in a matrix
 */
void clear_matrix(gf2matrix *m);
/**
 * \brief
 * sets a (row,col) bit to the given value
 */
void set_bit_at(gf2matrix *m, int bit, int row, int col);
/**
 * \brief
 * gets a bit locates at (row,col) in matrix m
 */
int get_bit_at(const gf2matrix *m, int row, int col);

/**
 * brief
 * extract a region of a given matrix and return it as another matrix
 * if dest is NULL then create a new matrix
 */
gf2matrix *extract_region(gf2matrix *dest, const gf2matrix *src,
		int start_row, int start_col,
		int row_count, int col_count);
/**
 * \brief
 * compare two regions in two matrices
 */
int comp_regions(const gf2matrix *a,
		int start_row_a, int start_col_a,
		const gf2matrix *b,
		int start_row_b, int start_col_b,
		int row_count, int col_count);

/**
 * \brief
 * populates a given matrix from bits of a given integer
 * until the count is reached or the matrix is fully populated
 */
int scalar2matrix(gf2matrix *dest, int scalar, int count);
/**
 * \brief
 * populates a given matrix from bits of a given integer
 * The scalar is copied from its least signifcant bit starting
 * from (offset_row, offser_col) + bitcount working
 * its way up until the offset is reached (row,col)
 */
int scalar2matrix_offset(gf2matrix *m, int scalar, int bitcount,
		int offset_row, int offset_col);
/**
 * \brief
 * populate a scalar from a given matrix starting at (0,0)
 * and up to bitcount number of elements
 */
int matrix2scalar(int *scalar, const gf2matrix *m, int bitcount);
/**
 * \brief
 * populate a scalar from a given matrix starting at offset
 * and up to bitcount number of elements
 */
int matrix2scalar_offset(int *scalar, const gf2matrix *m, int bitcount,
		int offset_row, int offset_col);
/**
 * copy src matrix to given top-left corrdinates in dest matrix
 */
int copy_matrix_to_offset(gf2matrix *dest, const gf2matrix *src,
		int dest_offset_row, int dest_offset_col);
/**
 * copy src matrix to the given destination matrix
 * src and dest must have identical dimensions
 */
int copy_matrix(gf2matrix *dest, const gf2matrix *src);
/**
 * \brief
 * print matrix to stdout
 */
void print_matrix(const gf2matrix *m, const char *label);
/**
 * \brief
 * populate a matrix with random values
 */
void randomize_matrix(gf2matrix *m);
/**
 * \brief
 * Calculate the rank of the given matrix
 */
int calc_rank(gf2matrix *m);

#endif /* _MATRIX_H_ */
