#include <stdint.h>
#include "m4ri/m4ri.h"
#include "m4ri_gf2matrix_typedef.h"
/*
 * File: m4ri_matrix_impl
 * Implementation of remaining matrix operations based on
 * M4RI library
 */
int copy_matrix_to_offset(gf2matrix *dest, const gf2matrix *src,
		int dest_offset_row, int dest_offset_col)
{
	size_t i, j;
	if (!src)
		return -1;
	for (i = 0; i < src->nrows; ++i) {
		for (j = 0; j < src->ncols; ++j) {
			BIT bit = mzd_read_bit(src, i, j);
			mzd_write_bit(dest, dest_offset_row + i, dest_offset_col + j, bit);
		}
	}
	return 0;
}
gf2matrix *new_matrix(int rows, int cols)
{
	return mzd_init(rows, cols);
}
void free_matrix(gf2matrix *m)
{
	if(!m)
		return;
	mzd_free(m);
}
gf2matrix *make_identity_matrix(int bits)
{
	mzd_t *I = mzd_init(bits, bits);
	if (I)
		mzd_set_ui(I, 1);
	return I;
}
void convert_to_identity(gf2matrix *m)
{
	if (!m)
		return;
	mzd_set_ui(m, 1);
}
gf2matrix *add_matrices(gf2matrix *C, const gf2matrix *A, const gf2matrix *B)
{
	if (!A || !B)
		return NULL;
	return mzd_add(C, (mzd_t *) A, (mzd_t *) B);
}
gf2matrix *mul_matrices(gf2matrix *C, const gf2matrix *A, const gf2matrix *B)
{
	if (!A || !B)
		return NULL;

	return mzd_mul_m4rm(C, (mzd_t *) A, (mzd_t *) B, 0);
}
int comp_matrices(const gf2matrix *a, const gf2matrix *b)
{
	if (!a || !b)
		return 1;
	return mzd_cmp(a, b);
}
int get_rows(const gf2matrix *m)
{
	if (!m)
		return 0;
	return m->nrows;
}
int get_cols(const gf2matrix *m)
{
	if (!m)
		return 0;
	return m->ncols;
}
gf2matrix *extract_region(gf2matrix *dest, const gf2matrix *src, int start_row,
		int start_col, int row_count, int col_count)
{
	if (!src)
		return NULL;
	return mzd_submatrix(dest, src, start_row, start_col,
			start_row + row_count, start_col + col_count);
}
int scalar2matrix_offset(gf2matrix *m, int scalar, int bitcount, int offset_row,
		int offset_col)
{
	int i, j;
	if (!m)
		return -1;
	/*
	 * start from the least significant bit and place
	 * it at the offset; then work your
	 * way up by going left and up starting
	 * from the offset column
	 */
	for (i = offset_row; i < m->nrows && bitcount > 0; ++i) {
		for (j = offset_col; j < m->ncols && bitcount > 0; ++j) {
			mzd_write_bit(m, i, j, (scalar >> (bitcount - 1)) & 0x1);
			--bitcount;
		}
	}
	return 0;
}
int scalar2matrix(gf2matrix *m, int scalar, int bitcount)
{
    return scalar2matrix_offset(m, scalar, bitcount, 0, 0);
}
int matrix2scalar_offset(int *scalar, const gf2matrix *m, int bitcount,
		int offset_row, int offset_col)
{
	int i, j;
	int rows = m->nrows;
	int cols = m->ncols;
	if (!m)
		return -1;
	if (bitcount > ((rows - offset_row) * (cols - offset_col)))
		return -1;
	/*
	 * start from the least significant bit and place
	 * it at the offset; then work your
	 * way up by going left and up starting
	 * from the offset column
	 */
	*scalar = 0;
	for (i = offset_row; i < rows && bitcount > 0; ++i) {
		for (j = offset_col; j < cols && bitcount > 0; ++j) {
			*scalar |= (mzd_read_bit(m, i, j) & 0x1) << (bitcount - 1);
			--bitcount;
		}
	}
	return 0;
}
int matrix2scalar(int *scalar, const gf2matrix *m, int bitcount)
{
    return matrix2scalar_offset(scalar, m, bitcount, 0, 0);
}
gf2matrix *dup_matrix(const gf2matrix *m)
{
	if (!m)
		return NULL;
	return mzd_copy(NULL, m);
}
void clear_matrix(gf2matrix *m)
{
	if (!m)
		return;
	mzd_set_ui(m, 0);
}
void set_bit_at(gf2matrix *m, int bit, int row, int col)
{
	if (!m)
		return;
	mzd_write_bit(m, row, col, bit);
}
int get_bit_at(const gf2matrix *m, int row, int col)
{
	if (!m)
		return 0;
	return (int) mzd_read_bit(m, row, col);
}
void print_matrix(const gf2matrix *m, const char *label)
{
	if (!m) {
		printf("%s\n[NULL]\n", label);
		return;
	}
	if (label)
		printf("%s\n", label);
	mzd_print(m);
}
void randomize_matrix(gf2matrix *m)
{
	mzd_randomize(m);
}
int calc_rank(const gf2matrix *m)
{
	gf2matrix *copy = dup_matrix(m);
	int r = (int) mzd_echelonize_m4ri(copy, 0, 0);
	free_matrix(copy);
	return r;
}
gf2matrix *invert_matrix(gf2matrix *dest, const gf2matrix *m)
{
    if (!m)
        return NULL;
    mzd_t *temp = mzd_inv_m4ri(NULL, (mzd_t*) m, 0);
    if(!temp) {
        dest = NULL;
    }
    else {
        mzd_t *I = make_identity_matrix(m->nrows);
        mzd_t *prod = mul_matrices(NULL, m, temp);
        if(comp_matrices(prod, I) == 0) {
            // inversion was successful
            if (dest) {
                mzd_copy(dest, temp);
                mzd_free(temp);
            }
            else {
                dest = temp;
            }
        }
        else {
            dest = NULL;
        }
        mzd_free(prod);
        mzd_free(I);
    }
    return dest;
}
