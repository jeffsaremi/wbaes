#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <stdint.h>
#include "nonlinear_ops.h"
#include "test.h"
#include "util.h"

extern void random_shuffle(uint8_t *p, int n);

static void print_sbox(const uint8_t *f, const int n)
{
	const int _2_n = 1 << n;
	int i, j;
	printf("f[] =   [");
	for (j = 0; j < _2_n; ++j)
		printf("%d ", f[j]);
	printf("]\n");
	for (i = n; i != 0; --i) {
		uint8_t fk[_2_n];
		extract_fk(fk, f, n, i);
		printf("fk[%d] = [", i);
		for (j = 0; j < _2_n; ++j)
			printf("%d ", fk[j]);
		printf("]\n");
	}
}

int test_is_sbox_sac()
{
	uint8_t f[] = { 3, 1, 4, 0, 2, 5, 6, 7 };
	ASSERT(is_sbox_sac(f, 3));
	print_sbox(f, 3);
}
int test_invert_sbox()
{
	int i, j, bits;
	for (bits = 3; bits < 8; ++bits) {
		int size = 1 << bits;
		uint8_t *sbox = malloc(size);
		uint8_t *invsbox = malloc(size);
		for (j = 0; j < size; ++j)
			sbox[j] = j;
		for (i = 0; i < 10; ++i) {
			random_shuffle(sbox, size);
			ASSERT(invert_sbox(invsbox, sbox, bits) == 0);
			/* print_sbox(invsbox,3); */
		}
		free(invsbox);
		free(sbox);
	}
}
int test_extend_sbox()
{
	uint8_t f[] = { 3, 1, 4, 0, 2, 5, 6, 7 };
	uint8_t f_nplus1[sizeof(f) * 2];
	uint8_t g[] = { 1, 0, 0, 0, 1, 1, 0, 1 };
	uint8_t fg[sizeof(f) * 2];
	unsigned int i, j, k;
	uint8_t expected[16] = { 11, 1, 4, 0, 10, 13, 6, 15, 9, 3, 8, 12, 5, 2, 7,
			14 };

	ASSERT(is_bit_vector_sac(g, 3));
	ASSERT(extend_f_with_g(fg, g, f, 3, 1));
	printf("fg[] =   [");
	for (i = 0; i < sizeof(f_nplus1); ++i)
		printf("%d ", fg[i]);
	printf("]\n");

	ASSERT(memcmp(fg, expected, sizeof(fg)) == 0);
	ASSERT(is_sbox_sac(fg, 4));

	for (k = 1; k < 4; ++k) {
		extend_f_by1(f_nplus1, f, 3, k);
		printf("f_nplus1[] =   [");
		for (j = 0; j < sizeof(f_nplus1); ++j)
			printf("%d ", f_nplus1[j]);
		printf("]\n");
		ASSERT(is_sbox_sac(f_nplus1, 4));
	}
}
int test_make_random_sbox()
{
	int i;
	for (i = 4; i < 8; ++i) {
		const int size = 1 << i;
		uint8_t *sbox = malloc(size);
		uint8_t *invsbox = malloc(size);
		ASSERT(make_random_sbox(sbox, invsbox, i) == 0);
		dump_hex("sbox ", sbox, size);
		free(sbox);
		free(invsbox);
	}
}

int main()
{
	srand(1234);

	TEST(test_is_sbox_sac);
	TEST(test_invert_sbox);
	TEST(test_extend_sbox);
	TEST(test_make_random_sbox);
	return 0;
}
