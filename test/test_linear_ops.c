#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include "test.h"
#include "linear_ops.h"

/*
 * checks to see if the given matrix is:
 * a) invertible as a whole, and
 * b) invertible block by block
 */
static int is_block_invertible(const gf2matrix *m)
{
	int rc = 0;
	int rows = get_rows(m);
	int cols = get_cols(m);
	int blocksize = 4;
	int i, j;
	gf2matrix *sub = NULL;
	if (!is_invertible(m))
		goto cleanup;
	/* return early if rows or cols not multiples of 2 or not square matrix */
	if ((rows % 4 != 0) || (cols % 4 != 0) || rows != cols)
		return 0;
	for (i = 0; i < rows && !rc; i += blocksize) {
		for (j = 0; j < cols && !rc; j += blocksize) {
			sub = extract_region(NULL, m, i, j, blocksize, blocksize);
			if (!sub) {
				printf("incorrect region specs: rows:%d, cols:%d,"
					" blocksize:%d, i:%d, j:%d\n", rows, cols, blocksize, i, j);
				goto cleanup;
			}
			if (!is_invertible(sub)) {
				print_matrix(sub, "non-invertible block: ");
				goto cleanup;
			}
			free_matrix(sub);
			sub = NULL;
		}
	}
	rc = 1;
	cleanup: if (sub)
		free_matrix(sub);
	return rc;
}
static int is_block_identity(const gf2matrix *m)
{
	int rc = 0;
	int rows = get_rows(m);
	int cols = get_cols(m);
	int blocksize = 4, i, j;
	gf2matrix *sub = NULL;
	gf2matrix *I = make_identity_matrix(4);
	/* return early if rows or cols not multiples of 2 or not square matrix */
	if ((rows % 2 != 0) || (cols % 2 != 0) || rows != cols)
		return 0;
	for (i = 0; i < rows && !rc; i += blocksize) {
		for (j = 0; j < cols && !rc; j += blocksize) {
			sub = extract_region(NULL, m, i, j, blocksize, blocksize);
			if (!sub) {
				printf("incorrect region specs: rows:%d, cols:%d,"
					" blocksize:%d, i:%d, j:%d\n", rows, cols, blocksize, i, j);
				goto cleanup;
			}
			if (comp_matrices(I, sub) != 0) {
				print_matrix(sub, "non-identity block: ");
				goto cleanup;
			}
			free_matrix(sub);
			sub = NULL;
		}
	}
	rc = 1;
	cleanup: if (sub)
		free_matrix(sub);
	free_matrix(I);
	return rc;
}
void test_free_matrices()
{
    int i;
    int n = 3;
    gf2matrix **m = (gf2matrix **) malloc(n * sizeof(gf2matrix*));
    for (i = 0; i < n; ++i) {
        m[i] = new_matrix(3, 3);
    }
    free_matrices(m, n);
    free(m);
}
void test_byte_vector_ops()
{
	int i, offset;
	gf2matrix *v = NULL;
	uint8_t byte;
	v = new_vector(16);
	for (offset = 0; offset < 8; ++offset) {
		for (i = 0; i < 256; ++i) {
			clear_matrix(v);
			ASSERT(byte2vector_offset(v, i, offset) == 0);
			ASSERT(vector2byte_offset(&byte, v, offset) == 0);
			ASSERT_EQUAL(i, byte);
		}
	}
	free_vector(v);
}
void test_matrix_mul_128x128()
{
	gf2matrix *m = new_matrix(128, 128);
	uint8_t bytes[16];
	uint8_t result[16], result2[16];
	int i;

	randomize_matrix(m);
	for (i = 0; i < 16; ++i) {
		bytes[i] = rand();
	}
	mul_array_by_matrix_128x128(result, m, bytes);
	/*  comparative section */
	{
		gf2matrix *v = new_vector(128);
		gf2matrix *p = NULL;
		for (i = 0; i < 16; ++i) {
			byte2vector_offset(v, bytes[i], i * 8);
		}
		p = mul_matrices(NULL, m, v);
		for (i = 0; i < 16; ++i) {
			vector2byte_offset(&result2[i], p, i * 8);
		}
		free_vector(v);
		free_matrix(p);
	}
	ASSERT(memcmp(result, result2, 16) == 0);
	free_matrix(m);
}
void test_matrix_mul_8x128()
{
	gf2matrix *m = new_matrix(128, 8);
	uint8_t byte = rand();
	uint8_t result[16], result2[16];
	int i;

	randomize_matrix(m);
	mul_byte_by_matrix_128x8(result, m, byte);
	/*  comparative section */
	{
		gf2matrix *v = new_vector(8);
		gf2matrix *p = NULL;

		byte2vector(v, byte);
		p = mul_matrices(NULL, m, v);
		for (i = 0; i < 16; ++i) {
			vector2byte_offset(&result2[i], p, i * 8);
		}
		free_vector(v);
		free_matrix(p);
	}
	ASSERT(memcmp(result, result2, 16) == 0);
	free_matrix(m);
}
void test_matrix_mul_32x32()
{
	gf2matrix *m = new_matrix(32, 32);
	uint8_t bytes[4];
	uint8_t result[4], result2[4];
	int i;

	randomize_matrix(m);
	for (i = 0; i < 4; ++i) {
		bytes[i] = rand();
	}
	mul_array_by_matrix_32x32(result, m, bytes);
	/*  comparative section */
	{
		gf2matrix *v = new_vector(32);
		gf2matrix *p = NULL;
		for (i = 0; i < 4; ++i) {
			byte2vector_offset(v, bytes[i], i * 8);
		}
		p = mul_matrices(NULL, m, v);
		for (i = 0; i < 4; ++i) {
			vector2byte_offset(&result2[i], p, i * 8);
		}
		free_vector(v);
		free_matrix(p);
	}
	ASSERT(memcmp(result, result2, 4) == 0);
	free_matrix(m);
}
void test_matrix_mul()
{
	gf2matrix *m = new_matrix(8, 8);
	uint8_t byte = rand();
	uint8_t result, result2;

	randomize_matrix(m);
	mul_byte_by_matrix(&result, m, byte);
	/*  comparative section */
	{
		gf2matrix *v = new_vector(8);
		gf2matrix *p = NULL;

		byte2vector(v, byte);
		p = mul_matrices(NULL, m, v);
		vector2byte(&result2, p);
		free_vector(v);
		free_matrix(p);
	}
	ASSERT_EQUAL(result , result2);
	free_matrix(m);
}
void test_extend_block_invertible()
{
	int bits[] = { 4, 8, 32, 128, -1 };
	int i;
	for (i = 0; bits[i] != -1; ++i) {
		gf2matrix *m, *minv;
		int rc = make_block_invertible_matrix_pair(&m, &minv, bits[i]);
		print_matrix(m, "m:");
		ASSERT(rc==0);
		/* ASSERT(is_block_invertible(m)); */
		free_matrix(m);
		free_matrix(minv);
	}
}
void test_extend_block_invertible2()
{
	int bits = 32;
	int i;
	gf2matrix *m, *minv;
	int rc = make_block_invertible_matrix_pair(&m, &minv, bits);
	print_matrix(m, "m:");
	ASSERT(rc==0);
	/* ASSERT(is_block_invertible(m)); */
	{
		uint8_t bytes[4], result[4], result2[4];
		for (i = 0; i < 4; ++i) {
			bytes[i] = rand();
		}
		mul_array_by_matrix_32x32(result, m, bytes);
		dump_hex("bytes: ", bytes, 4);
		dump_hex("result: ", result, 4);
		mul_array_by_matrix_32x32(result2, minv, result);
		dump_hex("result2: ", result2, 4);
	}
	free_matrix(m);
	free_matrix(minv);
}
void test_sliced_mul()
{
	gf2matrix *m;
	gf2matrix *slices[16];
	int i, j, k;
	uint8_t v[16][16];
	uint8_t input[16], result[16], expected[16];

	make_block_invertible_matrix_pair(&m, NULL, 128);
	slice_matrix_vertically(slices, 16, m);

	for (i = 0; i < 16; ++i) {
		input[i] = rand();
		mul_byte_by_matrix_128x8(v[i],
				slices[i], input[i]);
	}
	for(i=0; i < 16; ++i) {
		result[i] = v[15][i];
		for(k=14; k != -1; --k)
			result[i] = v[k][i] ^ result[i];
	}
	mul_array_by_matrix_128x128(expected, m, input);
/* 	dump_hex("input:\t", input, 16); */
	dump_hex("result:\t", result, 16);
	dump_hex("expected:\t", expected, 16);
	ASSERT(memcmp(expected, result, 16) == 0);
	free_matrices(slices, 16);
	free_matrix(m);
}
void test_sliced_mul_with_inv()
{
	gf2matrix *m, *minv;
	gf2matrix *slices[16], *slices_inv[16];
	int i, j, k;
	uint8_t v[16][16];
	uint8_t input[16], result[16], expected[16];

	make_block_invertible_matrix_pair(&m, &minv, 128);
	slice_matrix_vertically(slices, 16, m);
	slice_matrix_vertically(slices_inv, 16, minv);

	for (i = 0; i < 16; ++i) {
		input[i] = rand();
		mul_byte_by_matrix_128x8(v[i],
				slices[i], input[i]);
	}
	for(i=0; i < 16; ++i) {
		result[i] = v[15][i];
		for(k=14; k != -1; --k)
			result[i] = v[k][i] ^ result[i];
	}
	/*  now the inverse calculations */
	for (i = 0; i < 16; ++i) {
		mul_byte_by_matrix_128x8(v[i],
				slices_inv[i], result[i]);
	}
	for(i=0; i < 16; ++i) {
		expected[i] = v[15][i];
		for(k=14; k != -1; --k)
			expected[i] = v[k][i] ^ expected[i];
	}
	dump_hex("result:\t", input, 16);
	dump_hex("expected:\t", expected, 16);
	ASSERT(memcmp(expected, input, 16) == 0);
	free_matrices(slices, 16);
	free_matrices(slices_inv, 16);
	free_matrix(m);
	free_matrix(minv);
}
int main()
{
	/* srand(1234); */
    TEST(test_free_matrices);
	TEST(test_byte_vector_ops);
	TEST(test_matrix_mul_128x128);
	TEST(test_matrix_mul_8x128);
	TEST(test_matrix_mul_32x32);
	TEST(test_matrix_mul);
	TEST(test_extend_block_invertible);
	TEST(test_sliced_mul);
	TEST(test_sliced_mul_with_inv);
	return 0;
}
