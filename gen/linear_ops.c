#include <stdio.h>
#include <stdlib.h>
#include "linear_ops.h"
#include "util.h"

#define MAX_2X2_INVERTIBLES 6
#define MAX_A2_MATRICES 3

static gf2matrix *g_all_2x2_invertibles[MAX_2X2_INVERTIBLES][2];
static int created2x2invertibles = 0;

/**
 * \brief
 * populate a matrix made of sub-matrices I (identity) and 0 (zero)
 * r is the size of the identity sub-matrix
 * |1 0 0 0|
 * |0 1 0 0|
 * |0 0 0 0|
 * |0 0 0 0|
 * In the above example r = 2
 */
static void init_IO_matrix(gf2matrix *m, int r)
{
	int i;
	clear_matrix(m);
	for (i = 0; i < r; ++i)
		set_bit_at(m, 1, i, i);
}
/**
 * \brief
 * initialize and populate all 2 by 2 invertible matrices
 * A lot of other functions depend on this
 */
static void init_all_2x2_invertibles()
{
	/*
	 * |1 0|  |1 1|  |0 1|  |0 1|  |1 1|  |1 0|
	 * |0 1|  |1 0|  |1 1|  |1 0|  |0 1|  |1 1|
	 */
	const int all_invertible_scalars[MAX_2X2_INVERTIBLES][2] = {
			{ 0x9, 0x9 } /* identity */, { 0xE, 0x7 } /* 1110 */,
			{ 0x7, 0xE } /* 0111 */, { 0x6, 0x6 } /* 0110 */,
			{ 0xB, 0xB } /* 1011 */, { 0xD, 0xD } /* 1101 */, };
	if (!created2x2invertibles) {
		int i;
		for (i = 0; i < MAX_2X2_INVERTIBLES; ++i) {
			g_all_2x2_invertibles[i][0] = new_matrix(2, 2);
			g_all_2x2_invertibles[i][1] = new_matrix(2, 2);
			scalar2matrix(g_all_2x2_invertibles[i][0],
					all_invertible_scalars[i][0], 4);
			scalar2matrix(g_all_2x2_invertibles[i][1],
					all_invertible_scalars[i][1], 4);
		}
		for (i = 0; i < MAX_2X2_INVERTIBLES; ++i) {
/* 			printf("i: %d  ", i); */
/* 			print_matrix(g_all_2x2_invertibles[i][0], "2x2:"); */
		}
		created2x2invertibles = 1;
	}
}
static void get_2x2invertible_pair(int index, gf2matrix **m, gf2matrix **minv)
{
	init_all_2x2_invertibles();
	assert(index < MAX_2X2_INVERTIBLES && index >= 0);
	/* printf("index = %d\n", index); */
	*m = g_all_2x2_invertibles[index][0];
	assert(*m);
	if (minv)
		*minv = g_all_2x2_invertibles[index][1];
}
/**
 * \brief
 * return the A(bitsxbits) matrix corresponding to the given rank
 * see Generating Large Non-Singular ...; Jmes Xiao, section 2.
 * @param p number of bits
 * @param r rank
 * @return a matrix A of dims (pxp)
 *
 * For a 2x2 we have:
 *  r=0    r=1    r=2
 * |1 0|  |1 1|  |0 1|
 * |0 1|  |1 0|  |1 1|
 *
 * In general:
 * 1. r = 0 A = I(p,p)
 * 2. r is even:
 * 		 	    | |0 1|                        |
 *              | |1 1|                        |
 *              |      ...(r/2 times           |
 *      A =     |          diagonally)         |
 *              |                              |
 *              |                   I(p-r,p-r) |
 *
 * 3. r is 1:
 * 		 	    | |1 1|           |
 *     A =      | |1 0|           |
 *              |      I(p-2,p-2) |
 *
 * 4. r is odd and > 1:
 *              | |1 1 1|                      |
 *              | |1 1 0|                      |
 *              | |1 0 0|                      |
 * 		 	    |       |0 1|                  |
 *     A =      |       |1 1|                  |
 *              |            ...(n times       |
 *              |            diagonally)       |
 *              |                              |
 *              |                 I(p-2n,p-2n) |
 *    where n = (r - 3) / 2
 */
static gf2matrix *make_A_matrix(int p, int r)
{
	gf2matrix *a;
	if (r == 0) {
		a = make_identity_matrix(p);
	}
	else if (r == 1) {
		gf2matrix *I = make_identity_matrix(p - 2);
		gf2matrix *block = new_matrix(2, 2);
		scalar2matrix(block, 0xE, 4);
		a = new_matrix(p, p);
		clear_matrix(a);
		copy_matrix_to_offset(a, I, 2, 2);
		copy_matrix_to_offset(a, block, 0, 0);
		free_matrix(block);
		free_matrix(I);
	}
	else if ((r % 2) == 0) {
		int i;
		gf2matrix *I = NULL;
		gf2matrix *block = new_matrix(2, 2);
		if (p > r)
			I = make_identity_matrix(p - r);
		scalar2matrix(block, 0xE, 4);
		a = new_matrix(p, p);
		clear_matrix(a);
		if (I)
			copy_matrix_to_offset(a, I, r, r);
		for (i = 0; i < r; i += 2)
			copy_matrix_to_offset(a, block, i, i);
		free_matrix(block);
		free_matrix(I);
	}
	else /*  r is odd and > 1 */
	{
		int i;
		int n;
		gf2matrix *I = NULL;
		gf2matrix *block = NULL;
		n = (r - 3) / 2;
		if (n = 0) {
			a = make_identity_matrix(p);
		}
		else {
			int p_2n = p - (2 * n);
			block = new_matrix(2, 2);
			if (p > r)
				I = make_identity_matrix(p - r);
			scalar2matrix(block, 0xE, 4);
			a = new_matrix(p, p);
			clear_matrix(a);
			if (I)
				copy_matrix_to_offset(a, I, r, r);
			for (i = 0; i < n; ++i) {
				int offset = (i * 2) + 3;
				copy_matrix_to_offset(a, block, offset, offset);
			}
			set_bit_at(a, 1, 0, 0);
			set_bit_at(a, 1, 0, 1);
			set_bit_at(a, 1, 0, 2);
			set_bit_at(a, 1, 1, 1);
			set_bit_at(a, 1, 1, 0);
			set_bit_at(a, 1, 2, 0);
		}
		free_matrix(block);
		free_matrix(I);
	}
	/* print_matrix(a, "A:"); */
	return a;
}
static gf2matrix *make_A2_matrix(int r)
{
	gf2matrix *a = new_matrix(2, 2);
	const int A2_scalars[3] = { 0x9 /* identity */, 0xE /* 1110 */, 0x7 /* 0111 */
	};
	scalar2matrix(a, A2_scalars[r], 4);
	/* print_matrix(a, "A:"); */
	return a;
}
/**
 * \brief
 * Extract and copy a 2 by 2 block starting at the given start row
 * from the given src
 */
static gf2matrix *extract_block_row(gf2matrix *block, int start_row,
		const gf2matrix *src)
{
	return extract_region(block, src, start_row, 0, 2, get_cols(src));
}
/**
 * \brief
 * Extract and copy a 2 by 2 block starting at the given start column
 * from the given src
 */
static gf2matrix *extract_block_col(gf2matrix *block, int start_col,
		const gf2matrix *src)
{
	return extract_region(block, src, 0, start_col, get_rows(src), 2);
}
/**
 * \brief combine A, B, C, D into a new matrix
 * in the following order:
 *  | A  B |
 *  | C  D |
 */
static gf2matrix *combine_matrices(const gf2matrix *A, const gf2matrix *B,
		const gf2matrix *C, const gf2matrix *D)
{
	gf2matrix *result = NULL;
	assert(get_rows(A) + get_rows(C) == get_rows(B) + get_rows(D));
	assert(get_cols(A) + get_cols(B) == get_cols(C) + get_cols(D));
	result = new_matrix(get_rows(A) + get_rows(C), get_cols(A) + get_cols(B));
	assert(copy_matrix_to_offset(result, A, 0, 0) ==0);
	assert(copy_matrix_to_offset(result, C, get_rows(A), 0) ==0);
	assert(copy_matrix_to_offset(result, B, 0, get_cols(A)) ==0);
	assert(copy_matrix_to_offset(result, D, get_rows(A), get_cols(A)) ==0);
	return result;
}
/*
 **
 * \brief
 * Create the following from the given components:
 * |  M           Y             |
 * |  X  X.Minv.Y + Pinv.A.Qinv |
 */
static gf2matrix *extend_by_blocks(const gf2matrix *m, const gf2matrix *x,
		const gf2matrix *y, const gf2matrix *xminvy, const gf2matrix *pinv,
		const gf2matrix *a, const gf2matrix *qinv)
{
	gf2matrix *result = NULL;
	gf2matrix *addterm = NULL;
	/* print_matrix(xminvy, "X.Minv.Y: "); */
	{
		gf2matrix *aqinv = NULL, *pinvaqinv = NULL;
		aqinv = mul_matrices(NULL, a, qinv);
		assert(aqinv);
		pinvaqinv = mul_matrices(NULL, pinv, aqinv);
		assert(pinvaqinv);
		/* print_matrix(pinvaqinv, "Pinv.A.Qinv: "); */
		addterm = add_matrices(NULL, xminvy, pinvaqinv);
		free_matrix(aqinv);
		free_matrix(pinvaqinv);
	}
	assert(addterm);
	/* print_matrix(addterm, "X.Minv.Y + Pinv.A.Qinv: "); */
	assert(is_invertible(addterm));
	result = combine_matrices(m, y, x, addterm);
	assert(result);

	if (!is_invertible(result)) {
		/* printf("extended matrix was not invertible. incorrect calculations\n"); */
		/* print_matrix(result, "result: "); */
		free_matrix(result);
		result = NULL;
	}
	free_matrix(addterm);
	return result;
}
/**
 * @detail
 * All variables are named as in the fore mentioned paper
 * 1. Use a row of blocks of Matrix M to create a matrix X
 * 2. Use a column of blocks of Matrix M to create a matrix Y
 * 3. Get invertible matrices P and Q such that P(X.Minv.Y)Q = [I 0]
 *    (I is of r dimension; 0 is of 2-r dimension; 0 < r < 3)
 * 4. Define matrix A as follows:
 * 		a) if r = 0 then A = I
 * 		b) if r = 1 then A = [[1 1] [1 0]]
 * 		c) if r = 2 then A = [[0 1] [1 1]]
 * 5. Then the following is a (t+2, 2) block invertible matrix:
 * 		|  M           Y             |
 * 		|  X  X.Minv.Y + Pinv.A.Qinv |
 */
static gf2matrix *extend_block_invertible_by2(const gf2matrix *m)
{
	gf2matrix *result = NULL;
	gf2matrix *x = NULL, *y = NULL;
	int multiples = get_rows(m) / 2;
	int x_start = get_random(0, multiples - 1) * 2;
	int y_start = get_random(0, multiples - 1) * 2;
	if (!m) {
		gf2matrix *_2x2;
		get_2x2invertible_pair(get_random(0, MAX_2X2_INVERTIBLES - 1), &_2x2,
				NULL);
		result = dup_matrix(_2x2);
		return result;
	}
	/* print_matrix(m, "Extending m by 2:"); */

	/*  step 1 */
	x = extract_block_row(NULL, x_start, m);
	assert(x);
	/* print_matrix(x, "X:"); */

	/*  step 2 */
	y = extract_block_col(NULL, y_start, m);
	assert(y);
	/* print_matrix(y, "Y:"); */
	{ /*  steps 3, 4, 5 */
		int r; /* r is the rank of X.Minv.Y */
		gf2matrix *i0mat = new_matrix(2, 2);
		gf2matrix *prod = new_matrix(2, 2);
		gf2matrix *temp = NULL;
		gf2matrix *mInv = invert_matrix(NULL, m);
		gf2matrix *xminvy = NULL;
		gf2matrix *a2 = NULL;
		{
			assert(mInv);
			/*  calculate X.Minv.Y */
			gf2matrix *minvy = mul_matrices(NULL, mInv, y);
			xminvy = mul_matrices(NULL, x, minvy);
			assert(xminvy);
			free_matrix(minvy);
		}
		/* print_matrix(xminvy, "xM-1y: "); */
		r = calc_rank(xminvy);
		/* printf("rank of X.Minv.Y = %d\n", r); */
		init_IO_matrix(i0mat, r);
		/* print_matrix(i0mat, "I0:"); */
		a2 = make_A2_matrix(r);

		temp = new_matrix(2, 2);
		while (!result) {
			int i, j, i_tries, j_tries;
			for (i = get_random(0, MAX_2X2_INVERTIBLES - 1), i_tries = 0; i_tries
					< MAX_2X2_INVERTIBLES && !result; i = (i + 1)
					% MAX_2X2_INVERTIBLES, ++i_tries) {
				gf2matrix *q, *qinv;
				get_2x2invertible_pair(i, &q, &qinv);
				/* print_matrix(q, "Q:"); */
				mul_matrices(temp, xminvy, q);
				assert(temp);
				for (j = get_random(0, MAX_2X2_INVERTIBLES - 1), j_tries = 0; j_tries
						< MAX_2X2_INVERTIBLES && !result; j = (j + 1)
						% MAX_2X2_INVERTIBLES, ++j_tries) {
					gf2matrix *p, *pinv;
					get_2x2invertible_pair(j, &p, &pinv);
					mul_matrices(prod, p, temp);
					assert(prod);
					/* print_matrix(p, "P:"); */
					/* print_matrix(prod, "P.X.Minv.Y.Q:"); */
					if (comp_matrices(prod, i0mat) == 0) {
						/*  step 5 */
						result = extend_by_blocks(m, x, y, xminvy, pinv, a2,
								qinv);
						assert(result);
						break;
					}
				}
			}
		}
		free_matrix(a2);
		free_matrix(xminvy);
		free_matrix(mInv);
		free_matrix(prod);
		free_matrix(temp);
		free_matrix(i0mat);

		if (!result)
			printf("incorrect matrix expansion algorithm");
	}
	cleanup: free_matrix(y);
	free_matrix(x);
	return result;
}
/**
 * \brief
 * populate a matrix with a tile matrix
 */
static void tile_matrix(gf2matrix *out, const gf2matrix *tile)
{
	int c = get_cols(tile);
	int r = get_rows(tile);
	int i, j;
	for (i = 0; i < get_rows(out); i += r) {
		for (j = 0; j < get_cols(out); j += c) {
			copy_matrix_to_offset(out, tile, i, j);
		}
	}
}

int make_block_invertible_matrix_pair(gf2matrix **m, gf2matrix **minv, int bits)
{
	int rc = 0;
	int i, j;
	gf2matrix *out = NULL, *outinv = NULL;
	assert(bits % 4 == 0);
	assert(m);

	while (out == NULL || get_rows(out) != bits) {
		gf2matrix *temp = extend_block_invertible_by2(out);
		free_matrix(out);
		out = dup_matrix(temp);
		free_matrix(temp);
	}
	*m = out;
	if (minv) {
		outinv = invert_matrix(NULL, out);
		assert(outinv);
		*minv = outinv;
	}
	return rc;
}
int mul_byte_by_matrix(uint8_t *result, const gf2matrix *m, uint8_t byte)
{
	/**
	 * 1. convert byte to a matrix (vector)
	 * 2. multiply vector by m; obtain result matrix
	 * 3. convert result matrix into byte array
	 */
	int i;
	int octets = get_rows(m) / 8;
	gf2matrix *vec = new_matrix(8, 1);
	byte2vector(vec, byte);
	gf2matrix *prod = mul_matrices(NULL, (gf2matrix *) m, vec);
	for (i = 0; i < octets; ++i) {
		vector2byte_offset(&result[i], prod, i * 8);
	}
	free_matrix(vec);
	free_matrix(prod);
	return 0;
}
int mul_array_by_matrix(uint8_t *result, const gf2matrix *m,
		const uint8_t *bytes)
{
	/**
	 * 1. convert the input byte array to a matrix (vector)
	 * 2. multiply vector by m; obtain result matrix
	 * 3. convert result matrix into byte array
	 */
	int i;
	int octets = get_rows(m) / 8;
	gf2matrix *prod;
	gf2matrix *vec = new_matrix(get_rows(m), 1);
	for (i = 0; i < octets; ++i) {
		byte2vector_offset(vec, bytes[i], i * 8);
	}
	prod = mul_matrices(NULL, m, vec);
	for (i = 0; i < octets; ++i) {
		vector2byte_offset(&result[i], prod, i * 8);
	}
	free_matrix(vec);
	free_matrix(prod);
	return 0;
}
int mul_byte_by_matrix_128x8(uint8_t result[16], const gf2matrix *m,
		uint8_t byte)
{
	return mul_byte_by_matrix(result, m, byte);
}
int mul_array_by_matrix_128x128(uint8_t result[16], const gf2matrix *m,
		const uint8_t value[16])
{
	return mul_array_by_matrix(result, m, &value[0]);
}
int mul_byte_by_matrix_32x8(uint8_t *result, const gf2matrix *m, uint8_t byte)
{
	return mul_byte_by_matrix(result, m, byte);
}
int mul_array_by_matrix_32x32(uint8_t result[4], const gf2matrix *m,
		const uint8_t value[4])
{
	return mul_array_by_matrix(result, m, &value[0]);
}
int slice_matrix_vertically(gf2matrix *slices[], int count, const gf2matrix *src)
{
	int i;
	const size_t slice_cols = get_cols(src) / count;
	for (i = 0; i < count; ++i) {
		size_t start_col = i * slice_cols;
		slices[i] = extract_region(NULL, src, 0, start_col, get_rows(src),
				slice_cols);
	}
	return 0;
}
gf2matrix **new_sliced_matrix_vertical(int count, const gf2matrix *src)
{
	int i;
	const size_t slice_cols = get_cols(src) / count;
	gf2matrix **slices = (gf2matrix **) malloc(count * sizeof(gf2matrix *));
	for (i = 0; i < count; ++i) {
		slices[i] = new_matrix(get_rows(src), slice_cols);
	}
	slice_matrix_vertically(slices, count, src);
	return slices;
}
void free_matrices(gf2matrix **matrices, int count)
{
	while (count > 0) {
		free_matrix(matrices[count - 1]);
		--count;
	}
}
gf2matrix *new_vector(int rows)
{
	return new_matrix(rows, 1);
}
void free_vector(gf2matrix *v)
{
	free_matrix(v);
}
int byte2vector_offset(gf2matrix *v, uint8_t byte, int offset_row)
{
	if (!v)
		return -1;
	return scalar2matrix_offset(v, byte, 8, offset_row, 0);
}
int byte2vector(gf2matrix *v, uint8_t byte)
{
	return byte2vector_offset(v, byte, 0);
}
int vector2byte_offset(uint8_t *byte, const gf2matrix *v, int offset_row)
{
	int a, rc;
	if (!v)
		return -1;
	rc = matrix2scalar_offset(&a, v, 8, offset_row, 0);
	*byte = (uint8_t) a;
	return rc;
}
int vector2byte(uint8_t *byte, const gf2matrix *v)
{
	return vector2byte_offset(byte, v, 0);
}
