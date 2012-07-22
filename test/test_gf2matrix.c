#include <stdio.h>
#include <stdlib.h>
#include "gf2matrix.h"
#include "util.h"
#include "test.h"

#ifndef MIN
#define MIN(a,b) ((a) < (b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#endif

void test_matrix_allocation()
{
	gf2matrix *m = new_matrix(3, 3);
	free_matrix(m);
}
void test_get_rows_cols()
{
	int rows = 3;
	int cols = 4;
	gf2matrix *m = new_matrix(rows, cols);
	ASSERT(rows == get_rows(m));
	ASSERT(cols == get_cols(m));
	free_matrix(m);
}
void test_identity_matrix()
{
	int rows = 2;
	gf2matrix *m = make_identity_matrix(rows);
	ASSERT(1 == get_bit_at(m, 0, 0));
	ASSERT(0 == get_bit_at(m, 0, 1));
	ASSERT(0 == get_bit_at(m, 1, 0));
	ASSERT(1 == get_bit_at(m, 1, 1));

	free_matrix(m);
}
void test_convert_to_identity()
{
	int rows = 2;
	int cols = 2;
	gf2matrix *m = new_matrix(rows, cols);
	convert_to_identity(m);
	ASSERT(1 == get_bit_at(m, 0, 0));
	ASSERT(0 == get_bit_at(m, 0, 1));
	ASSERT(0 == get_bit_at(m, 1, 0));
	ASSERT(1 == get_bit_at(m, 1, 1));

	free_matrix(m);
}
void test_scalar2matrix()
{
	gf2matrix *m = new_matrix(2, 2);
	ASSERT(scalar2matrix(m, 0xE /* [[1 1] [1 0]] */, 4) == 0);
	/* print_matrix(m, "m: "); */
	ASSERT(1 == get_bit_at(m, 0, 0));
	ASSERT(1 == get_bit_at(m, 0, 1));
	ASSERT(1 == get_bit_at(m, 1, 0));
	ASSERT(0 == get_bit_at(m, 1, 1));
	free_matrix(m);
}
void test_scalar2matrix_offset()
{
	int row = 2;
	int col = 4;
	int count = 4;
	gf2matrix *m = new_matrix(4, 6);
	int scalar = get_random(1, 0xf);
	/* printf("scalar = %x\n", scalar); */
	ASSERT(scalar2matrix_offset(m, scalar, count,
					row, col) == 0);
	/* print_matrix(m, "m: "); */
	ASSERT(((scalar >> (count-1)) & 0x1) == get_bit_at(m, row, col));
	--count;
	ASSERT(((scalar >> (count-1)) & 0x1) == get_bit_at(m, row, col+1));
	--count;
	ASSERT(((scalar >> (count-1)) & 0x1) == get_bit_at(m, row+1, col));
	--count;
	ASSERT(((scalar >> (count-1)) & 0x1) == get_bit_at(m, row+1, col+1));
	free_matrix(m);
}
void test_matrix2scalar()
{
	int n, scalar;
	int i;
	for (i = 0; i < 10; ++i) {
		gf2matrix *m = NULL;
		int bits = get_random(2, 16);
		int max_n = 0;
		int j;
		for (j = 0; j < bits - 1; ++j) {
			max_n |= (1 << j);
		}
		m = new_matrix((bits + 1) / 2, (bits + 1) / 2);
		/* printf("rows=%d, cols=%d\n", get_rows(m), get_cols(m)); */
		clear_matrix(m);
		n = get_random(1, max_n);
		ASSERT(scalar2matrix(m, n, bits) == 0);
		/* print_matrix(m, "m:"); */
		/* printf("n=0x%x, bits = %d\n", n, bits); */
		ASSERT(matrix2scalar(&scalar, m, bits) == 0);
		/* printf("scalar=0x%x\n", scalar); */
		ASSERT(scalar == n);
		free_matrix(m);
	}
}
void test_matrix2scalar_offset()
{
	int n, scalar;
	int i;
	for (i = 0; i < 10; ++i) {
		gf2matrix *m = NULL;
		int bits = get_random(2, 16);
		int offset;
		int max_n = 0;
		int j;
		for (j = 0; j < bits - 1; ++j) {
			max_n |= (1 << j);
		}
		m = new_matrix(bits * 2, bits * 2);
		offset = get_random(0, bits);
		/* printf("rows=%d, cols=%d\n", get_rows(m), get_cols(m)); */
		clear_matrix(m);
		n = get_random(1, max_n);
		ASSERT(scalar2matrix_offset(m, n, bits, offset, offset) == 0);
		/* print_matrix(m, "m:"); */
		/* printf("n=0x%x, bits = %d\n", n, bits); */
		ASSERT(matrix2scalar_offset(&scalar, m, bits, offset, offset) == 0);
		/* printf("scalar=0x%x\n", scalar); */
		ASSERT(scalar == n);
		free_matrix(m);
	}
}
void test_is_invertible()
{
	int rows = 2;
	gf2matrix *m = new_matrix(rows, rows);
	gf2matrix *I = make_identity_matrix(rows);
	scalar2matrix(m, 0xA /* [[1 0] [1 0]] */, 4);
	ASSERT(1 == is_invertible(I));
	ASSERT(0 == is_invertible(m));

	free_matrix(I);
	free_matrix(m);
}
void test_invert()
{
	int rows = 2;
	gf2matrix *m = NULL, *expected = NULL, *minv = NULL;
	m = new_matrix(rows, rows);
	expected = new_matrix(rows, rows);
	scalar2matrix(m, 0xE /* [[1 1] [1 0]] */, 4);
	scalar2matrix(expected, 0x7 /* [[0 1] [1 1]] */, 4);
	minv = invert_matrix(NULL, m);
	ASSERT(minv);
	ASSERT(comp_matrices(minv, expected) == 0);
	free_matrix(minv);
	free_matrix(expected);
	free_matrix(m);
}
void test_dup()
{
	int rows = 2;
	gf2matrix *m = NULL, *dupped = NULL;
	m = new_matrix(rows, rows);
	scalar2matrix(m, 0xE /* [[1 1] [1 0]] */, 4);
	dupped = dup_matrix(m);
	ASSERT(dupped);
	ASSERT(comp_matrices(m, dupped) == 0);
	free_matrix(dupped);
	free_matrix(m);
}
void test_add()
{
	int rows = 2;
	gf2matrix *a = NULL, *b = NULL, *aplusb = NULL, *added = NULL;
	a = new_matrix(rows, rows);
	b = new_matrix(rows, rows);
	aplusb = new_matrix(rows, rows);
	scalar2matrix(a, 0xE /* [[1 1] [1 0]] */, 4);
	scalar2matrix(b, 0x7 /* [[0 1] [1 1]] */, 4);
	scalar2matrix(aplusb, 0x9 /* [[1 0] [0 1]] */, 4);
	added = add_matrices(NULL, a, b);
	/*
	 print_matrix(aplusb, "aplusb:");
	 print_matrix(added, "added:");
	 */
	ASSERT(added);
	ASSERT(comp_matrices(aplusb, added) == 0);
	free_matrix(added);
	free_matrix(aplusb);
	free_matrix(b);
	free_matrix(a);
}
void test_mul()
{
	int rows = 2;
	gf2matrix *a = NULL, *b = NULL, *ab = NULL, *mul = NULL;
	a = new_matrix(rows, rows);
	b = new_matrix(rows, rows);
	ab = new_matrix(rows, rows);
	scalar2matrix(a, 0xE /* [[1 1] [1 0]] */, 4);
	scalar2matrix(b, 0x7 /* [[0 1] [1 1]] */, 4);
	scalar2matrix(ab, 0x9 /* [[1 0] [0 1]] */, 4);
	mul = mul_matrices(NULL, a, b);
	/*
	 print_matrix(ab, "ab:");
	 print_matrix(mul, "mul:");
	 */
	ASSERT(mul);
	ASSERT(comp_matrices(ab, mul) == 0);
	free_matrix(mul);
	free_matrix(ab);
	free_matrix(b);
	free_matrix(a);
}
void test_clear()
{
	int rows = 3;
	int i, j;
	gf2matrix *m = make_identity_matrix(rows);
	clear_matrix(m);
	for (i = 0; i < rows; ++i) {
		for (j = 0; j < rows; ++j) {
			ASSERT(0 == get_bit_at(m, i, j));
		}
	}
	free_matrix(m);
}
void test_get_set_bit()
{
	int rows = 2;
	int cols = 4;
	int i, j;
	gf2matrix *m = new_matrix(rows, cols);
	for (i = 0; i < rows; ++i) {
		for (j = 0; j < cols; ++j) {
			set_bit_at(m, 1, i, j);
			ASSERT(1 == get_bit_at(m, i, j));
		}
	}
	free_matrix(m);
}
void test_extract_region()
{
	int rows = 4;
	int cols = 6;
	int i, j;
	int row_count, col_count;
	gf2matrix *m = new_matrix(rows, cols);
	randomize_matrix(m);
	/* print_matrix(m, "m:"); */
	for (i = rows - 1; i > 0; --i) {
		for (j = cols - 1; j > 0; --j) {
			row_count = MAX(1, rows-i);
			col_count = MAX(1, cols-j);
			gf2matrix *r = extract_region(NULL, m, i, j, row_count, col_count);
			/* print_matrix(r, "r:"); */
			ASSERT(0 == comp_regions(m, i, j,
							r, 0, 0,
							row_count, col_count));
			free_matrix(r);
		}
	}
	free_matrix(m);
}
void test_copy()
{
	int bigrows = 8;
	int smallrows = 2;
	int rowoffset, coloffset;
	gf2matrix *big = new_matrix(bigrows, bigrows);
	gf2matrix *small = new_matrix(smallrows, smallrows);
	randomize_matrix(big);
	randomize_matrix(small);

	for (rowoffset = 0; rowoffset < bigrows - smallrows; ++rowoffset) {
		for (coloffset = 0; coloffset < bigrows - smallrows; ++coloffset) {
			copy_matrix_to_offset(big, small, rowoffset, coloffset);
			ASSERT(0 == comp_regions(big, rowoffset, coloffset,
							small, 0, 0,
							smallrows, smallrows));
		}
	}
	free_matrix(small);
	free_matrix(big);
}

int main()
{
	TEST(test_matrix_allocation);
	TEST(test_get_rows_cols);
	TEST(test_identity_matrix);
	TEST(test_convert_to_identity);
	TEST(test_scalar2matrix);
	TEST(test_scalar2matrix_offset);
	TEST(test_matrix2scalar);
	TEST(test_matrix2scalar_offset);
	TEST(test_is_invertible);
	TEST(test_invert);
	TEST(test_dup);
	TEST(test_add);
	TEST(test_mul);
	TEST(test_clear);
	TEST(test_get_set_bit);
	TEST(test_extract_region);
	TEST(test_copy);
	return 0;
}
