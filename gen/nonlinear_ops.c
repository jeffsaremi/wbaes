#include <stdlib.h>
#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include "nonlinear_ops.h"
#include "util.h"

/**
 * \brief
 * randomly shuffles an array of chars of size n
 */
void random_shuffle(uint8_t *p, int n)
{
	uint8_t temp;
	int i;
	for (i = n - 1; i > 0; --i) {
		int r = get_random(0, i);
		temp = p[r];
		p[r] = p[i];
		p[i] = temp;
	}
}
int is_bit_vector_sac(const uint8_t *f, int n)
{
	const int _2_n = 1 << n;
	const int _2_n_1 = 1 << (n - 1);
	int x, k;
	for (k = 1; k <= n; ++k) {
		int Ck = 1 << (k - 1);
		int sum_Gf2_n = 0;
		for (x = 0; x < _2_n; ++x) {
			int temp = f[x] ^ f[x ^ Ck];
			sum_Gf2_n += temp;
			sum_Gf2_n %= _2_n;
		} /*  loop for x in Z2n */
		/* 		printf("expected=%d, sum=%d, 2^n-1=%d, n=%d, k=%d\n", */
		/* 				_2_n_1 * _2_n, sum_Gf2_n, _2_n_1, n, k); */
		if (sum_Gf2_n != _2_n_1)
			return 0;
	} /*  loop for Ci */
	return 1;
}
int is_sbox_sac(const uint8_t *f, int n)
{
	const int _2_n = 1 << n;
	const int _2_n_1 = 1 << (n - 1);
	int i, x;
	for (i = 1; i <= n; ++i) {
		int Ci = 1 << (i - 1);
		int sum_Gf2_n = 0;
		for (x = 0; x < _2_n; ++x) {
			int temp = f[x] ^ f[x ^ Ci];
			sum_Gf2_n += temp;
			sum_Gf2_n %= _2_n;
		} /*  loop for x in Z2n */
		/* printf("sum=%d, 2^n-1=%d\n", sum_Gf2_n, _2_n_1); */
		if (sum_Gf2_n != _2_n_1)
			return 0;
	} /*  loop for Ci */
	return 1;
}
uint8_t *extend_f_with_g(uint8_t *dest, const uint8_t *g, const uint8_t *f,
		int n, int k)
{
	const int _2_n = 1 << n;
	const int _2_nplus1 = _2_n << 1;
	const int Ck = 1 << (k - 1);
	int D1_g[_2_nplus1];
	int x;
	for (x = 0; x < _2_n; ++x) {
		D1_g[x] = g[x]; /*  D1_g */
		dest[x] = f[x]; /*  D0_f */
	}
	for (x = _2_n; x < _2_nplus1; ++x) {
		uint8_t temp = (x - _2_n) ^ Ck;
		D1_g[x] = D1_g[temp] ^ 1; /*  D1_g */
		dest[x] = f[temp]; /*  D0_f */
	}
	for (x = 0; x < _2_nplus1; ++x) {
		/* printf("D1_g[%d] = %d, dest[%d]=%d\n ", x, D1_g[x], x, dest[x]); */
		dest[x] = (dest[x] | (D1_g[x] << n)); /*  % _2_nplus1; */
	}
	return dest;
}
uint8_t *extract_fk(uint8_t *dest, const uint8_t *f, int n, int k)
{
	const int _2_n = 1 << n;
	const int Ck = 1 << (k - 1);
	int x;
	for (x = 0; x < _2_n; ++x) {
		dest[x] = (f[x] & Ck) == 1;
	}
	return dest;
}
uint8_t *extend_f_by1(uint8_t *dest, const uint8_t *f, int n, int k)
{
	/*  get f_k */
	/*  make Dk1_f_k */
	/*  make Dk0_f */
	/*  concat the above */
	const int _2_n = 1 << n;
	uint8_t D1_fk[_2_n];

	/* printf("2^n=%d, 2^n+1=%d\n", _2_n, _2_nplus1); */
	extract_fk(D1_fk, f, n, k);
	return extend_f_with_g(dest, D1_fk, f, n, k);
}
int invert_sbox(uint8_t *invsbox, const uint8_t *sbox, int n)
{
	int rc = 0;
	const int _2_n = 1 << n;
	int x;
	assert(n <= 8); /*  only uint8_t size */
	assert(invsbox);

	for (x = 0; x < _2_n; ++x) {
		invsbox[sbox[x]] = x;
		/* 		printf("x=%d, sbox[%d]=%d, invsbox[%d]=%d\n", */
		/* 				x, x, sbox[x], sbox[x], invsbox[sbox[x]]); */
	}
	/* printf("n = %d\n", n); */
	/* dump_decimal("inverted sbox: ", invsbox, _2_n); */
	/* dump_decimal("sbox: ", sbox, _2_n); */
	for (x = 0; x < _2_n; ++x) {
		if (x != invsbox[sbox[x]]) {
			printf("non-invertible\n");
			rc = -1;
			break;
		}
	}
	return rc;
}
int make_random_sbox(uint8_t *sbox, uint8_t *invsbox, int bits)
{
	int rc = 0;
	uint8_t _3bitSbox[8] = { 0, 1, 2, 3, 4, 5, 6, 7 };
	const int _2_n = 1 << bits;
	uint8_t temp1[_2_n];
	uint8_t temp2[_2_n];
	int i;
	assert(bits >= 3);
	assert(sbox);
	assert(invsbox);

	while (1) {
		random_shuffle(_3bitSbox, 8);
		/* dump_decimal("shuffled ", _3bitSbox, 8); */
		if (is_sbox_sac(_3bitSbox, 3)) {
			uint8_t *src;
			uint8_t *dest;
			memcpy(temp1, _3bitSbox, 8);
			for (i = 3; i < bits; ++i) {
				src = &temp1[0];
				dest = &temp2[0];
				int k = get_random(1, i);
				/* printf("k=%d\n", k); */
				/* dump_decimal("pre extension ", src, (1 << i)); */
				extend_f_by1(dest, src, i, k);
				/* dump_decimal("post extension ", dest, (1 << (i+1))); */
				assert(is_sbox_sac(dest, i+1));

				memcpy(temp1, temp2, (1 << (i + 1)));
			}
			memcpy(sbox, temp2, _2_n);
			assert(invert_sbox(invsbox, sbox, bits) == 0);
			/* 			if(!sbox || !invsbox || */
			/* 					invert_sbox(invsbox, sbox, bits) != 0) */
			/* 				rc = -1; */
			break;
		}
	}
	return rc;
}

