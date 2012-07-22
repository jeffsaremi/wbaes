#include <stdio.h>
#include <memory.h>
#include <assert.h>
#include "wb_aes_gen.h"
#include "linear_ops.h"
#include "nonlinear_ops.h"
#include "aes.h"
#include "util.h"

typedef int (*shift_col_f)(int, int);
typedef void (*mix_slices_f)(uint8_t *, uint8_t, int);

static const int DIRECTION_FORWARD = 0;
static const int DIRECTION_INVERSE = 1;

static uint8_t g_shifted_cols[4][4] = {
		{ 0, 1, 2, 3 },
		{ 1, 2, 3, 0 },
		{ 2, 3, 0, 1 },
		{ 3, 0, 1, 2 } };
static uint8_t g_inv_shifted_cols[4][4] = {
		{ 0, 1, 2, 3 },
		{ 3, 0, 1, 2 },
		{ 2, 3, 0, 1 },
		{ 1, 2, 3, 0 } };

int get_shifted_col(int row, int col)
{
	return g_shifted_cols[row][col];
}
int get_inv_shifted_col(int row, int col)
{
	return g_inv_shifted_cols[row][col];
}

void sub_bytes_hi_lo32(uint8_t v[4], uint8_t sbox[][16])
{
	int i;
	for (i = 0; i < 4; ++i) {
		v[i] = sbox[((i * 2) + 1)][(v[i] & 0x0F)] |
				(sbox[i * 2][(v[i] >> 4)] << 4);
	}
}
void sub_bytes_hi_lo128(uint8_t v[16], uint8_t sbox[][16])
{
	int i;
	for (i = 0; i < 16; ++i) {
		v[i] = sbox[((i * 2) + 1)][(v[i] & 0x0F)] |
				(sbox[i * 2][(v[i] >> 4)] << 4);
	}
}
int make_tbox(tbox_t tbox,
		uint8_t sbox[256],
		uint32_t expanded_key[(NR+1)*4],
		tbox_mixing_bijections_t tbox_mixing_bijection)
{
	int round, col, row, x, k;
	int rc = 0;
	for (round = 0; round < NR - 1; ++round) {
		uint8_t *key = (uint8_t *) &expanded_key[(round) * 4];
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				uint8_t shifted_index = (uint8_t) (row + 4
						* get_shifted_col(row, col));
				for (x = 0; x < 256; ++x) {
					uint8_t xpresub;
					mul_byte_by_matrix(
							&xpresub,
							tbox_mixing_bijection[round][row][get_shifted_col(row, col)],
							x);
					tbox[round][row][col][x] = sbox[(xpresub
							^ key[shifted_index])];
				}
			}
		}
	}
	{ /*  The last two rounds are merged into one */
		round = NR - 1;
		uint8_t *key = (uint8_t *) &expanded_key[(round) * 4];
		uint8_t *lastKey = (uint8_t *) &expanded_key[NR * 4];
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				uint8_t shifted_index = (uint8_t) (row + 4
						* get_shifted_col(row, col));
				uint8_t lastKey_index = (uint8_t) (row + (col * 4));
				for (x = 0; x < 256; ++x) {
					uint8_t xpresub;
					mul_byte_by_matrix(
							&xpresub,
							tbox_mixing_bijection[round][row][get_shifted_col(row, col)],
							x);
					tbox[round][row][col][x] = sbox[(xpresub
							^ key[shifted_index])] ^ lastKey[lastKey_index];
				}
			}
		}
	}
	return rc;
}
int make_inv_tbox(tbox_t tbox,
		uint8_t sbox[256],
		uint32_t expanded_key[(NR+1)*4],
		tbox_mixing_bijections_t tbox_mixing_bijection)
{
	int round, col, row, x, k;
	int rc = 0;
	for (round = NR - 1; round != 0; --round) {
		uint8_t *key = (uint8_t *) &expanded_key[(round + 1) * 4];
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				uint8_t shifted_index = (uint8_t) (row + 4
						* get_inv_shifted_col(row, col));
				for (x = 0; x < 256; ++x) {
					uint8_t xpresub;
					mul_byte_by_matrix(
							&xpresub,
							tbox_mixing_bijection[round][row][get_inv_shifted_col(row, col)],
							x);
					tbox[round][row][col][x] = sbox[(xpresub
							^ key[shifted_index])];
				}
			}
		}
	}
	{ /*  The last two rounds are merged into one */
		round = 0;
		uint8_t *key = (uint8_t *) &expanded_key[(round + 1) * 4];
		uint8_t *lastKey = (uint8_t *) &expanded_key[0];
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				uint8_t shifted_index = (uint8_t) (row + 4
						* get_inv_shifted_col(row, col));
				uint8_t lastKey_index = (uint8_t) (row + (col * 4));
				for (x = 0; x < 256; ++x) {
					uint8_t xpresub;
					mul_byte_by_matrix(
							&xpresub,
							tbox_mixing_bijection[round][row][get_inv_shifted_col(row, col)],
							x);
					tbox[round][row][col][x] = sbox[(xpresub
							^ key[shifted_index])] ^ lastKey[lastKey_index];
				}
			}
		}
	}
	return rc;
}
int make_typeII(typeII_t typeII,
		tbox_t tbox,
		gf2matrix *mix_columns_mixing_bijection,
		sboxes_8bit_t decoding_sbox[NR-1],
		sboxes_32bit_t encoding_sbox[NR-1])
{
	int round, col, row, x, k;
	int rc = 0;
	/*  there's a ShiftRows step before TypeII step */
	/*  hence the decoding sboxes need to be adjusted */
	for (round = 0; round < (NR - 1); ++round) {
		for (col = 0; col < 4; ++col) {
			for (row = 0; row < 4; ++row) {
				for (x = 0; x < 256; ++x) {
					uint8_t slice[4] = {0,0,0,0};
					uint8_t xmapped = sub_bytes_hi_lo( x,
							(decoding_sbox[round])[row][get_shifted_col(row, col)]);
					calc_mixing_slices(slice, tbox[round][row][col][xmapped], row);
					rc = mul_array_by_matrix_32x32(slice,
							mix_columns_mixing_bijection, slice);
					assert(rc == 0);
					if(rc != 0)
						goto end;
					for (k = 0; k < 4; ++k) {
						/*  care must be taken to use the proper encoding sbox */
						/*  for results; Note that items from one slice are */
						/*  not being XOR'd but rather corresponding items */
						/*  from slices generated for one column */
						slice[k] = sub_bytes_hi_lo( slice[k],
							(uint8_t (*)[16])&(encoding_sbox[round])[k][col][2 * row]);
					}
					memcpy(typeII[round][row][col][x], slice, 4);
				}
			}
		}
	}
end:
	return rc;
}
int make_inv_typeII(typeII_t typeII,
		tbox_t tbox,
		gf2matrix *mix_columns_mixing_bijection,
		sboxes_8bit_t decoding_sbox[NR-1],
		sboxes_32bit_t encoding_sbox[NR-1])
{
	int round, col, row, x, k;
	int rc = 0;
	/*  there's a ShiftRows step before TypeII step */
	/*  hence the decoding sboxes need to be adjusted */
	for (round = NR - 2; round != -1; --round) {
		for (col = 0; col < 4; ++col) {
			for (row = 0; row < 4; ++row) {
				for (x = 0; x < 256; ++x) {
					uint8_t slice[4] = {0,0,0,0};
					uint8_t xmapped = sub_bytes_hi_lo( x,
							(decoding_sbox[round])[row][get_inv_shifted_col(row, col)]);
					calc_inv_mixing_slices(slice, tbox[round + 1][row][col][xmapped], row);
					rc = mul_array_by_matrix_32x32(slice,
							mix_columns_mixing_bijection, slice);
					assert(rc == 0);
					if(rc != 0)
						goto end;
					for (k = 0; k < 4; ++k) {
						/*  care must be taken to use the proper encoding sbox */
						/*  for results; Note that items from one slice are */
						/*  not being XOR'd but rather corresponding items */
						/*  from slices generated for one column */
						slice[k] = sub_bytes_hi_lo( slice[k],
							(uint8_t (*)[16])&(encoding_sbox[round])[k][col][2 * row]);
					}
					memcpy(typeII[round][row][col][x], slice, 4);
				}
			}
		}
	}
end:
	return rc;
}
int make_rounds_typeIV(typeIV32_t typeIV32,
		sboxes_32bit_t decoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox_inv[NR - 1])
{
	int round, col, row, x, y, pair, bitset;
	int rc = 0;
	for (round = 0; round < (NR - 1); ++round) {
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				for (pair = 0; pair < 3; ++pair) {
					for (bitset = 0; bitset < 2; ++bitset) {
						uint8_t *x_decoding_sbox =
								&(decoding_sbox)[round][row][col][pair * 2
										+ bitset][0];
						uint8_t *y_decoding_sbox =
								&(encoding_sbox_inv)[round][row][col][bitset][0];
						/*  only the last pair uses the decoding corresponding to the previous round */
						if (pair == 2) {
							y_decoding_sbox
									= &decoding_sbox[round][row][col][((pair + 1) * 2) + bitset][0];
						}
						for (x = 0; x < 16; ++x) {
							for (y = 0; y < 16; ++y) {
								uint8_t x_decoded = x_decoding_sbox[x];
								uint8_t y_decoded = y_decoding_sbox[y];
								uint8_t xy_xored = x_decoded ^ y_decoded;
								int xy_index = (x << 4) | y;
								if (xy_index > 127) {
									typeIV32[round][row][col][pair][bitset][xy_index % 128]
											|= (encoding_sbox[round][row][col][bitset][xy_xored] << 4);
								}
								else {
									typeIV32[round][row][col][pair][bitset][xy_index]
											= encoding_sbox[round][row][col][bitset][xy_xored];
								}
							}
						}
					}
				}
			}
		}
	}
	return rc;
}

int make_typeIV128(typeIV128_t typeIV128,
		sboxes_128bit_t decoding_sbox,
		sboxes_8bit_t encoding_sbox,
		sboxes_8bit_t encoding_sbox_inv)
{
	int col, row, x, y, pair, bitset;
	int rc = 0;
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			for (pair = 0; pair < 15; ++pair) {
				for (bitset = 0; bitset < 2; ++bitset) {
					uint8_t *x_decoding_sbox =
							&decoding_sbox[row][col][pair * 2 + bitset][0];
					uint8_t *y_decoding_sbox =
							&encoding_sbox_inv[row][col][bitset][0];
					/*  only the last pair uses the decoding corresponding to the previous round */
					if (pair == 14) {
						y_decoding_sbox = &decoding_sbox[row][col][((pair + 1) * 2) + bitset][0];
					}
					for (x = 0; x < 16; ++x) {
						for (y = 0; y < 16; ++y) {
							uint8_t x_decoded = x_decoding_sbox[x];
							uint8_t y_decoded = y_decoding_sbox[y];
							uint8_t xy_xored = x_decoded ^ y_decoded;
							int xy_index = (x << 4) | y;
							if (xy_index > 127) {
								typeIV128[row][col][pair][bitset][xy_index % 128]
										|= (encoding_sbox[row][col][bitset][xy_xored] << 4);
							}
							else {
								typeIV128[row][col][pair][bitset][xy_index]
										= encoding_sbox[row][col][bitset][xy_xored];
							}
						}
					}
				}
			}
		}
	}
	return rc;
}
int make_typeIV_IA(typeIV128_t typeIV_IA,
		sboxes_128bit_t decoding_sbox,
		sboxes_8bit_t encoding_sbox,
		sboxes_8bit_t encoding_sbox_inv)
{
	return make_typeIV128(typeIV_IA, decoding_sbox,
			encoding_sbox, encoding_sbox_inv);
}
int make_typeIV_IB(typeIV128_t typeIV_IB,
		sboxes_128bit_t decoding_sbox,
		sboxes_8bit_t encoding_sbox,
		sboxes_8bit_t encoding_sbox_inv)
{
	return make_typeIV128(typeIV_IB, decoding_sbox,
			encoding_sbox, encoding_sbox_inv);
}
int make_typeIA(typeIA_t typeIA,
		gf2matrix *first_inv_tbox_mixing_bijection[4][4],
		gf2matrix *initial_decoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox)
{
	int col, row, x, k;
	int rc = 0;
	gf2matrix *decoding_slices[4 * 4];
	slice_matrix_vertically(decoding_slices, 4 * 4, initial_decoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row * 4 + col;
			for (x = 0; x < 256; ++x) {
				uint8_t xmapped = sub_bytes_hi_lo(x, decoding_sbox[row][col]);
				uint8_t slice[16];
				mul_byte_by_matrix_128x8(slice, decoding_slices[index],
						xmapped);
				for (k = 0; k < 16; ++k) {
					mul_byte_by_matrix(&slice[k],
							first_inv_tbox_mixing_bijection[k/4][k%4],
							slice[k]);
				}
				for (k = 0; k < 16; ++k) {
					/*  care must be taken to use the proper encoding sbox */
					/*  for results; Note that items from one slice are */
					/*  not being XOR'd but rather corresponding items */
					/*  from slices generated for one column */
					slice[k] = sub_bytes_hi_lo(slice[k],
						(uint8_t (*)[16])&(encoding_sbox[k/4][k%4][2 * index]));
				}
				/* sub_bytes_hi_lo128(result, encoding_sbox[row][col]); */
				memcpy(typeIA[row][col][x], slice, 16);
			}
		}
	}
	free_matrices(decoding_slices, 4*4);
	return rc;
}
int make_typeIB_ex(typeIB_t typeIB,
		int direction,
		uint8_t last_round_tbox[4][4][256],
		gf2matrix *final_encoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox)
{
	/*  there is a ShiftRows step before the Type IB step */
	/*  decoding of the inputs should take that into account */
	int col, row, x, k;
	int rc = 0;
	shift_col_f shift_col = (direction == 0)
			? get_shifted_col
			: get_inv_shifted_col;
	gf2matrix *encoding_slices[4 * 4];
	slice_matrix_vertically(encoding_slices, 4 * 4, final_encoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row * 4 + col;
			for (x = 0; x < 256; ++x) {
				uint8_t xmapped = sub_bytes_hi_lo(x,
						decoding_sbox[row][shift_col(row, col)]);
				uint8_t slice[16];
				xmapped = last_round_tbox[row][col][xmapped];
				mul_byte_by_matrix_128x8(slice, encoding_slices[index],
						xmapped);
				for (k = 0; k < 16; ++k) {
					/*  care must be taken to use the proper encoding sbox */
					/*  for results; Note that items from one slice are */
					/*  not being XOR'd but rather corresponding items */
					/*  from slices generated for one column */
					slice[k] = sub_bytes_hi_lo(slice[k],
						(uint8_t (*)[16])&(encoding_sbox[k/4][k%4][2 * index]));
				}
/* 				sub_bytes_hi_lo128(result, encoding_sbox[row][col]); */
				memcpy(typeIB[row][col][x], slice, 16);
			}
		}
	}
	free_matrices(encoding_slices, 4*4);
	return rc;
}
int make_typeIB(typeIB_t typeIB,
		uint8_t last_round_tbox[4][4][256],
		gf2matrix *final_encoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox)
{
	return make_typeIB_ex(typeIB, DIRECTION_FORWARD, last_round_tbox,
			final_encoding, decoding_sbox, encoding_sbox);
}
int make_inv_typeIB(typeIB_t typeIB,
		uint8_t last_round_tbox[4][4][256],
		gf2matrix *final_encoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox)
{
	return make_typeIB_ex(typeIB, DIRECTION_INVERSE, last_round_tbox,
			final_encoding, decoding_sbox, encoding_sbox);
}
int make_typeIII(typeIII_t typeIII,
		gf2matrix *inv_mix_columns_mixing_bijection,
		tbox_mixing_bijections_t inv_tbox_mixing_bijections,
		sboxes_8bit_t decoding_sbox[NR-1],
		sboxes_32bit_t encoding_sbox[NR-1])
{
	int round, col, row, x, k;
	int rc = 0;
	for (round = 0; round < (NR - 1); ++round) {
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				for (x = 0; x < 256; ++x) {
					uint8_t xmapped = sub_bytes_hi_lo(x,
							decoding_sbox[round][row][col]);
					uint8_t mc[4] = { 0, 0, 0, 0 };
					mc[row] = xmapped;
					rc = mul_array_by_matrix_32x32(mc,
							inv_mix_columns_mixing_bijection,
							mc);
					assert(rc == 0);
					if(rc != 0)
						goto end;
					for (k = 0; k < 4; ++k) {
						rc = mul_byte_by_matrix(&mc[k],
								inv_tbox_mixing_bijections[round][k][col],
								mc[k]);
					}
					for (k = 0; k < 4; ++k) {
						/*  care must be taken to use the proper encoding sbox */
						/*  for results; Note that items from one slice are */
						/*  not being XOR'd but rather corresponding items */
						/*  from slices generated for one column */
						mc[k] = sub_bytes_hi_lo( mc[k],
								(uint8_t (*)[16]) encoding_sbox[round][k][col][2 * row]);
					}
					memcpy(typeIII[round][row][col][x], mc, 4);
				}
			}
		}
	}
end:
	return rc;
}
void calc_mixing_slices(uint8_t mc[4], uint8_t tbox_val, int slice_num)
{
	uint8_t v = tbox_val;
	uint8_t v2 = xtime(v);
	uint8_t v3 = v2 ^ v;
	switch (slice_num) {
	case 0:
		/*  slice 0 */
		mc[0] = v2;
		mc[1] = v;
		mc[2] = v;
		mc[3] = v3;
		break;
	case 1:
		/*  slice 1 */
		mc[0] = v3;
		mc[1] = v2;
		mc[2] = v;
		mc[3] = v;
		break;
	case 2:
		/*  slice 2 */
		mc[0] = v;
		mc[1] = v3;
		mc[2] = v2;
		mc[3] = v;
		break;
	case 3:
		/*  slice 3 */
		mc[0] = v;
		mc[1] = v;
		mc[2] = v3;
		mc[3] = v2;
		break;
	}
}
void calc_inv_mixing_slices(uint8_t mc[4], uint8_t tbox_val,
		int slice_num)
{
	uint8_t i1 = (uint8_t) tbox_val;
	uint8_t i2 = xtime(i1);
	uint8_t i4 = xtime(i2);
	uint8_t i8 = xtime(i4);
	uint8_t i9 = i8 ^ i1;
	uint8_t iB = i8 ^ i2 ^ i1;
	uint8_t iD = i8 ^ i4 ^ i1;
	uint8_t iE = i8 ^ i4 ^ i2;

	switch (slice_num) {
	case 0:
		/*  slice 0 */
		mc[0] = iE;
		mc[1] = i9;
		mc[2] = iD;
		mc[3] = iB;
		break;
	case 1:
		/*  slice 1 */
		mc[0] = iB;
		mc[1] = iE;
		mc[2] = i9;
		mc[3] = iD;
		break;
	case 2:
		/*  slice 2 */
		mc[0] = iD;
		mc[1] = iB;
		mc[2] = iE;
		mc[3] = i9;
		break;
	case 3:
		/*  slice 3 */
		mc[0] = i9;
		mc[1] = iD;
		mc[2] = iB;
		mc[3] = iE;
		break;
	}
}
int make_rounds_sbox_pair_32(sboxes_32bit_t encoding_sboxes[],
		sboxes_32bit_t decoding_sboxes[], int rounds)
{
	int round, row, col, i;
	int rc = 0;
	for (round = 0; round < rounds && !rc; ++round) {
		for (row = 0; row < 4 && !rc; ++row) {
			for (col = 0; col < 4 && !rc; ++col) {
				for (i = 0; i < 8 && !rc; ++i) {
					rc = make_random_sbox(
							encoding_sboxes[round][row][col][i],
							decoding_sboxes[round][row][col][i],
							4);
					assert(rc == 0);
				}
			}
		}
	}
	return rc;
}
int make_rounds_sbox_pair_8(sboxes_8bit_t encoding_sboxes[],
		sboxes_8bit_t decoding_sboxes[],
		int rounds)
{
	int round, row, col, i;
	int rc = 0;
	for (round = 0; round < rounds && !rc; ++round) {
		for (row = 0; row < 4 && !rc; ++row) {
			for (col = 0; col < 4 && !rc; ++col) {
				for (i = 0; i < 2 && !rc; ++i) {
					rc = make_random_sbox(
							encoding_sboxes[round][row][col][i],
							decoding_sboxes[round][row][col][i],
							4);
					assert(rc == 0);
				}
			}
		}
	}
	return rc;
}
int make_sbox_pair_128(sboxes_128bit_t encoding_sboxes,
		sboxes_128bit_t decoding_sboxes)
{
	int row, col, i;
	int rc = 0;
	for (row = 0; row < 4 && !rc; ++row) {
		for (col = 0; col < 4 && !rc; ++col) {
			for (i = 0; i < 32 && !rc; ++i) {
				rc = make_random_sbox(encoding_sboxes[row][col][i],
						decoding_sboxes[row][col][i], 4);
				assert(rc == 0);
			}
		}
	}
	return rc;
}
int make_sbox_pair_32(sboxes_32bit_t encoding_sboxes,
		sboxes_32bit_t decoding_sboxes)
{
	int row, col, i;
	int rc = 0;
	for (row = 0; row < 4 && !rc; ++row) {
		for (col = 0; col < 4 && !rc; ++col) {
			for (i = 0; i < 8 && !rc; ++i) {
				rc = make_random_sbox(encoding_sboxes[row][col][i],
						decoding_sboxes[row][col][i], 4);
				assert(rc == 0);
			}
		}
	}
	return rc;
}
int make_sbox_pair_8(sboxes_8bit_t encoding_sboxes,
		sboxes_8bit_t decoding_sboxes)
{
	int row, col, i;
	int rc = 0;
	for (row = 0; row < 4 && !rc; ++row) {
		for (col = 0; col < 4 && !rc; ++col) {
			for (i = 0; i < 2 && !rc; ++i) {
				rc = make_random_sbox(encoding_sboxes[row][col][i],
						decoding_sboxes[row][col][i], 4);
				assert(rc == 0);
			}
		}
	}
	return rc;
}
int make_typeIV_II(typeIV32_t typeIV_II,
		sboxes_32bit_t decoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox_inv[NR - 1])
{
	return make_rounds_typeIV(typeIV_II, decoding_sbox,
			encoding_sbox, encoding_sbox_inv);
}
int make_typeIV_III(typeIV32_t typeIV_III,
		sboxes_32bit_t decoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox_inv[NR - 1])
{
	return make_rounds_typeIV(typeIV_III, decoding_sbox,
			encoding_sbox, encoding_sbox_inv);
}
int make_tbox_mixing_bijections(tbox_mixing_bijections_t tbox_mbs,
		tbox_mixing_bijections_t tbox_mb_invs)
{
	int round, row, col;
	for(round = 0; round < NR; ++round) {
		for(row = 0; row < 4; ++row) {
			for(col = 0; col < 4; ++col) {
				if(!tbox_mb_invs)
					make_block_invertible_matrix_pair(
							&(tbox_mbs[round][row][col]),
							NULL,
							8);
				else
					make_block_invertible_matrix_pair(
							&(tbox_mbs[round][row][col]),
							&(tbox_mb_invs[round][row][col]),
							8);
			}
		}
	}
	return 0;
}
void free_tbox_mixing_bijections(tbox_mixing_bijections_t tbox_mbs)
{
	int round, row, col;
	for(round = 0; round < NR; ++round) {
		for(row = 0; row < 4; ++row) {
			for(col = 0; col < 4; ++col) {
				free_matrix(tbox_mbs[round][row][col]);
			}
		}
	}
}
