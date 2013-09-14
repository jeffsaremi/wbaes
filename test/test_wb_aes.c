#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdint.h>
#include "wb_aes.h"
#include "linear_ops.h"
#include "nonlinear_ops.h"
#include "wb_aes_gen.h"
#include "util.h"
#include "test.h"
#include "test_util.h"

extern uint8_t SBox[256];
extern uint8_t ISBox[256];
extern void shift_rows(uint8_t state[4][4]);

void test_xor_hi_lo()
{
	static const int hibits_lookup = 0;
	static const int lowbits_lookup = 1;
	uint8_t dual8to4_lookup[2][128];
	int x, y;
	for (x = 0; x < 8; ++x) {
		for (y = 0; y < 16; ++y) {
			dual8to4_lookup[hibits_lookup][(x << 4) | y] = x ^ y;
			dual8to4_lookup[lowbits_lookup][(x << 4) | y] = x ^ y;
		}
	}
	for (x = 8; x < 16; ++x) {
		for (y = 0; y < 16; ++y) {
			dual8to4_lookup[hibits_lookup][((x << 4) | y) % 128] |= ((x ^ y)
					<< 4);
			dual8to4_lookup[lowbits_lookup][((x << 4) | y) % 128] |= ((x ^ y)
					<< 4);
		}
	}
	for (x = 0; x < 256; ++x) {
		for (y = 0; y < 256; ++y) {
			ASSERT((x^y) == do_xor_hi_lo(x, y, dual8to4_lookup));
		}
	}
}
void test_mix_columns_mixing_bijection()
{
	gf2matrix *mix_columns_mixing_bijection, *inv_mix_columns_mixing_bijection;
	int round, col;
	uint8_t state[4][4];
	uint8_t state2[4][4];

	make_block_invertible_matrix_pair(&mix_columns_mixing_bijection,
			&inv_mix_columns_mixing_bijection, 32);
	for (round = 0; round < 10; ++round) {
		printf("verifying round %d\n", round);
		randomize_state(state);
		dump_state("State before ", state);
		for (col = 0; col < 4; ++col) {
			uint8_t temp[4];
			copy_state_col(temp, state, col);
			uint8_t slice0[4], slice1[4], slice2[4], slice3[4];
			calc_mixing_slices(slice0, temp[0], 0);
			calc_mixing_slices(slice1, temp[1], 1);
			calc_mixing_slices(slice2, temp[2], 2);
			calc_mixing_slices(slice3, temp[3], 3);
			mul_array_by_matrix_32x32(slice0, mix_columns_mixing_bijection,
					slice0);
			mul_array_by_matrix_32x32(slice1, mix_columns_mixing_bijection,
					slice1);
			mul_array_by_matrix_32x32(slice2, mix_columns_mixing_bijection,
					slice2);
			mul_array_by_matrix_32x32(slice3, mix_columns_mixing_bijection,
					slice3);

			temp[0] = slice0[0] ^ slice1[0] ^ slice2[0] ^ slice3[0];
			temp[1] = slice0[1] ^ slice1[1] ^ slice2[1] ^ slice3[1];
			temp[2] = slice0[2] ^ slice1[2] ^ slice2[2] ^ slice3[2];
			temp[3] = slice0[3] ^ slice1[3] ^ slice2[3] ^ slice3[3];

			copy_col2state(state2, temp, col);
		}
		dump_state("State after ", state2);
		/*  verify by doing an Inv MixColumn and multiplication by invTbox MB */
		for (col = 0; col < 4; ++col) {
			uint8_t temp[4];
			copy_state_col(temp, state2, col);
			mul_array_by_matrix_32x32(temp, inv_mix_columns_mixing_bijection,
					temp);

			uint8_t slice0[4], slice1[4], slice2[4], slice3[4];
			calc_inv_mixing_slices(slice0, temp[0], 0);
			calc_inv_mixing_slices(slice1, temp[1], 1);
			calc_inv_mixing_slices(slice2, temp[2], 2);
			calc_inv_mixing_slices(slice3, temp[3], 3);

			temp[0] = slice0[0] ^ slice1[0] ^ slice2[0] ^ slice3[0];
			temp[1] = slice0[1] ^ slice1[1] ^ slice2[1] ^ slice3[1];
			temp[2] = slice0[2] ^ slice1[2] ^ slice2[2] ^ slice3[2];
			temp[3] = slice0[3] ^ slice1[3] ^ slice2[3] ^ slice3[3];

			copy_col2state(state2, temp, col);
		}
		dump_state("State inverted ", state2);
		ASSERT(comp_states(state, state2));
	}
	free_matrix(mix_columns_mixing_bijection);
	free_matrix(inv_mix_columns_mixing_bijection);
}
void test_tbox()
{
	tbox_mixing_bijections_t tbox_mixing_bijections, inv_tbox_mixing_bijections;
	tbox_t tbox;
	uint32_t expanded_key[(NR + 1) * 4];
	uint8_t state[4][4];
	uint8_t state2[4][4];
	int round, row, col, i;

	make_tbox_mixing_bijections(tbox_mixing_bijections,
			inv_tbox_mixing_bijections);

	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand(); /*  identity key expansion */
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);

	for (round = 0; round < NR - 1; ++round) {
		printf("round: %d \n", round);
		randomize_state(state);
		memcpy(&state2[0][0], &state[0][0], 16); /*  for validation */

		dump_state("State before ", state);
		shift_rows(state);
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = tbox[round][row][col][(state[row][col])];
			}
		}
		dump_state("State after TBox ", state);
		/*  validation */
		/*  The above should be equal to: */
		/*   0. mix */
		/* 	1. addroundkey */
		/*   2. subbytes */
		/* 	3. shiftrows */

		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				mul_byte_by_matrix(&state2[row][col],
						tbox_mixing_bijections[round][row][col],
						state2[row][col]);
			}
		}
		add_round_key(state2, expanded_key, round);
		sub_bytes(state2, SBox);
		shift_rows(state2);

		dump_state("Validation State ", state2);
		ASSERT(comp_states(state, state2));
	}
	/*  validation for the last round is different */
	round = NR - 1;
	{
		printf("rounds 9 and 10\n");
		randomize_state(state);
		memcpy(&state2[0][0], &state[0][0], 16); /*  for validation */

		dump_state("State before ", state);
		shift_rows(state);
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = tbox[round][row][col][(state[row][col])];
			}
		}
		dump_state("State after TBox ", state);
		/*  validation */
		/*  The above should be equal to: */
		/*   0. mix */
		/* 	1. addroundkey */
		/*   2. subbytes */
		/* 	3. shiftrows */
		/* 	4. addroundkey */

		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				mul_byte_by_matrix(&state2[row][col],
						tbox_mixing_bijections[round][row][col],
						state2[row][col]);
			}
		}
		add_round_key(state2, expanded_key, round);
		sub_bytes(state2, SBox);
		shift_rows(state2);
		add_round_key(state2, expanded_key, NR); /*  the last key */

		dump_state("Validation State ", state2);
		ASSERT(comp_states(state, state2));
	}
	free_tbox_mixing_bijections(tbox_mixing_bijections);
	free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
}
void test_typeIA()
{
	gf2matrix *inv_tbox_mixing_bijection[4][4];
	int row, col, i, j, k;
	typeIA_t typeIAs;
	gf2matrix *initial_decoding;
	sboxes_8bit_t input_sbox, input_sbox_inv;
	sboxes_128bit_t output_sbox, output_sbox_inv;
	_4bit_strip128_t strips;
	_4bit_strip128_t strips2, strips_temp;
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	gf2matrix *decoding_slices[4 * 4];

	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			make_block_invertible_matrix_pair(
					&inv_tbox_mixing_bijection[row][col], NULL, 8);
		}
	}
	make_block_invertible_matrix_pair(&initial_decoding, NULL, 128);
	make_sbox_pair_8(input_sbox, input_sbox_inv);
	make_sbox_pair_128(output_sbox, output_sbox_inv);

	make_typeIA(typeIAs, inv_tbox_mixing_bijection, initial_decoding,
			input_sbox_inv, output_sbox);

	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			pre_enc_state[row][col] = rand();
			state[row][col] = sub_bytes_hi_lo(pre_enc_state[row][col],
					input_sbox[row][col]);
		}
	}
	/*  end of preparation */
	do_typeIA(strips, state, typeIAs);
	/*  additionally decode the output and present the state */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			sub_bytes_hi_lo128(strips[row][col],
					output_sbox_inv[row][col]);
		}
	}
	dump_4bit_strip128("output: ", strips);

	slice_matrix_vertically(decoding_slices, 4 * 4, initial_decoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			mul_byte_by_matrix_128x8(strips_temp[row][col],
					decoding_slices[row * 4 + col],
					pre_enc_state[row][col]);
			for (k = 0; k < 16; ++k) {
				mul_byte_by_matrix(&strips_temp[row][col][k],
						inv_tbox_mixing_bijection[k/4][k%4],
						strips_temp[row][col][k]);
			}
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			for(k = 0; k < 16; ++k)
				strips2[k/4][k%4][index] = strips_temp[row][col][k];
		}
	}
	free_matrices(decoding_slices, 4 * 4);
	free_matrix(initial_decoding);
	for (i = 0; i < 4; ++i)
		free_matrices(inv_tbox_mixing_bijection[i], 4);
	dump_4bit_strip128("comparative output: ", strips2);
	ASSERT(comp_4bit_strips128(strips, strips2));
}
void test_typeIB()
{
	tbox_mixing_bijections_t tbox_mixing_bijections;
	tbox_t tbox;
	typeIB_t typeIBs;
	sboxes_8bit_t input_sbox, input_sbox_inv;
	sboxes_128bit_t output_sbox, output_sbox_inv;
	uint32_t expanded_key[(NR + 1) * 4];
	_4bit_strip128_t strips;
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	gf2matrix *final_encoding;
	gf2matrix *encoding_slices[4 * 4];
	_4bit_strip128_t strips2, strips_temp;
	int row, col, i, k;

	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand();
	make_tbox_mixing_bijections(tbox_mixing_bijections, NULL);
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);
	make_block_invertible_matrix_pair(&final_encoding, NULL, 128);
	make_sbox_pair_8(input_sbox, input_sbox_inv);
	make_sbox_pair_128(output_sbox, output_sbox_inv);

	make_typeIB(typeIBs, tbox[NR - 1], final_encoding, input_sbox_inv,
			output_sbox);

	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			pre_enc_state[row][col] = rand();
			state[row][col] = sub_bytes_hi_lo(pre_enc_state[row][col],
					input_sbox[row][col]);
		}
	}
	/*  end of preparation */
	shift_rows(state);
	do_typeIB(strips, state, typeIBs);
	/*  additionally decode the output and present the state */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			sub_bytes_hi_lo128(strips[row][col], output_sbox_inv[row][col]);
		}
	}
	dump_4bit_strip128("strips: ", strips);

	/*  comparative section */
	shift_rows(pre_enc_state);
	slice_matrix_vertically(encoding_slices, 4 * 4, final_encoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			uint8_t tbox_val = tbox[NR - 1][row][col][pre_enc_state[row][col]];
			mul_byte_by_matrix_128x8(strips_temp[row][col], encoding_slices[row * 4
					+ col], tbox_val);
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			for(k = 0; k < 16; ++k)
				strips2[k/4][k%4][index] = strips_temp[row][col][k];
		}
	}
	dump_4bit_strip128("comparative output: ", strips2);
	ASSERT(comp_4bit_strips128(strips, strips2));
	free_matrices(encoding_slices, 4*4);
	free_matrix(final_encoding);
	free_tbox_mixing_bijections(tbox_mixing_bijections);
}
void test_typeIV128()
{
	typeIV128_t typeIV_Is;
	sboxes_128bit_t input_sbox, input_sbox_inv;
	sboxes_8bit_t output_sbox, output_sbox_inv;
	_4bit_strip128_t strips;
	uint8_t pre_enc_strips[4][4][16];
	uint8_t state[4][4];
	uint8_t state2[4][4];
	int row, col, i, j, k;

	make_sbox_pair_128(input_sbox, input_sbox_inv);
	make_sbox_pair_8(output_sbox, output_sbox_inv);

	make_typeIV128(typeIV_Is, input_sbox_inv, output_sbox, output_sbox_inv);
	/*  prepare encoded input to the step */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			for (k = 0; k < 16; ++k) {
				pre_enc_strips[i][j][k] = rand();
				strips[i][j][k] = sub_bytes_hi_lo(pre_enc_strips[i][j][k],
						(uint8_t(*)[16]) &input_sbox[i][j][k * 2][0]);
			}
		}
	}
	/*  end of preparation */
	do_typeIV128(state, strips, typeIV_Is);
	/*  additionally decode the output and present the state */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			state[row][col] = sub_bytes_hi_lo(state[row][col],
					output_sbox_inv[row][col]);
		}
	}
	dump_state("state: ", state);

	/*  equivalent transformation: xor the pre-encoded strips */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			state2[row][col] = pre_enc_strips[row][col][15];
			for (k = 14; k != -1; --k) {
				state2[row][col] = pre_enc_strips[row][col][k]
						^ state2[row][col];
			}
		}
	}
	dump_state("comparative state: ", state2);
	ASSERT(comp_states(state, state2));
}
void test_typeIV32()
{
	typeIV32_t typeIV_IIs;
	sboxes_32bit_t input_sbox[NR - 1], input_sbox_inv[NR - 1];
	sboxes_8bit_t output_sbox[NR - 1], output_sbox_inv[NR - 1];
	_4bit_strip32_t strips;
	_4bit_strip32_t pre_enc_strips;
	uint8_t state[4][4];
	uint8_t state2[4][4];
	int round, row, col, i, j;

	make_rounds_sbox_pair_32(input_sbox, input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_8(output_sbox, output_sbox_inv, NR - 1);
	make_typeIV_II(typeIV_IIs, input_sbox_inv, output_sbox, output_sbox_inv);

	for (round = 0; round < NR - 1; ++round) {
		printf("testing round %d\n", round);
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				for (j = 0; j < 4; ++j) {
					pre_enc_strips[row][col][j] = rand();
					strips[row][col][j] = sub_bytes_hi_lo(
									pre_enc_strips[row][col][j],
									(uint8_t(*)[16]) &input_sbox[round][row][col][j
											* 2][0]);
				}
			}
		}
		/*  end of preparation */
		do_typeIV_II(state, strips, typeIV_IIs[round]);
		/*  additionally decode the output and present the state */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = sub_bytes_hi_lo(state[row][col],
						output_sbox_inv[round][row][col]);
			}
		}
		dump_state("state: ", state);
		/*  comparative section */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state2[row][col] = pre_enc_strips[row][col][0]
						^ (pre_enc_strips[row][col][1]
								^ (pre_enc_strips[row][col][2]
										^ pre_enc_strips[row][col][3]));
			}
		}
		dump_state("comparative state: ", state2);
		ASSERT(comp_states(state, state2));
	}
}
void test_typeIA_IV_combined()
{
	gf2matrix *inv_tbox_mixing_bijection[4][4];
	typeIA_t typeIAs;
	typeIV_IA_t typeIV_Is;
	gf2matrix *initial_decoding;
	sboxes_8bit_t input_sbox, input_sbox_inv;
	sboxes_128bit_t middle_sbox, middle_sbox_inv;
	sboxes_8bit_t output_sbox, output_sbox_inv;
	_4bit_strip128_t strips;
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	uint8_t strips2[4][4][16];
	gf2matrix *decoding_slices[4 * 4];
	uint8_t state2[4][4];
	int row, col, i, j, k;

	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			make_block_invertible_matrix_pair(
					&inv_tbox_mixing_bijection[row][col], NULL, 8);
		}
	}
	make_block_invertible_matrix_pair(&initial_decoding, NULL, 128);
	make_sbox_pair_8(input_sbox, input_sbox_inv);
	make_sbox_pair_128(middle_sbox, middle_sbox_inv);
	make_sbox_pair_8(output_sbox, output_sbox_inv);

	make_typeIA(typeIAs, inv_tbox_mixing_bijection, initial_decoding,
			input_sbox_inv, middle_sbox);
	make_typeIV_IA(typeIV_Is, middle_sbox_inv, output_sbox, output_sbox_inv);

	/*  prepare encoded input to the step */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			pre_enc_state[i][j] = rand();
			state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
					input_sbox[i][j]);
		}
	}
	/* dump_state("intial state: ", state); */
	/* dump_state("pre_enc_state: ", pre_enc_state); */

	do_typeIA(strips, state, typeIAs);
	do_typeIV_IA(state, strips, typeIV_Is);
	/*  decode output and present */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			state[i][j] = sub_bytes_hi_lo(state[i][j], output_sbox_inv[i][j]);
		}
	}
	dump_state("state: ", state);

	/*  equivalent transformation */
	slice_matrix_vertically(decoding_slices, 4 * 4, initial_decoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			mul_byte_by_matrix_128x8(strips2[row][col],
					decoding_slices[row * 4 + col], pre_enc_state[row][col]);
			for (k = 0; k < 16; ++k) {
				mul_byte_by_matrix(&strips2[row][col][k],
						inv_tbox_mixing_bijection[k/4][k%4],
						strips2[row][col][k]);
			}
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			state2[row][col] = strips2[15/4][15%4][index];
			for (k = 14; k != -1; --k) {
				state2[row][col] = strips2[k/4][k%4][index] ^ state2[row][col];
			}
		}
	}
	dump_state("comparative state: ", state2);
	ASSERT(comp_states(state, state2));
	free_matrices(decoding_slices, 4 * 4);
	free_matrix(initial_decoding);
	for (i = 0; i < 4; ++i)
		free_matrices(inv_tbox_mixing_bijection[i], 4);
}
void test_typeIB_IV_combined()
{
	tbox_mixing_bijections_t tbox_mixing_bijections;
	typeIB_t typeIBs;
	typeIV_IB_t typeIV_Is;
	gf2matrix *final_encoding;
	tbox_t tbox;
	uint32_t expanded_key[(NR + 1) * 4];
	sboxes_8bit_t input_sbox, input_sbox_inv;
	sboxes_128bit_t middle_sbox, middle_sbox_inv;
	sboxes_8bit_t output_sbox, output_sbox_inv;
	_4bit_strip128_t strips;
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	uint8_t state2[4][4];
	gf2matrix *encoding_slices[4 * 4];
	int row, col, i, j, k;

	make_tbox_mixing_bijections(tbox_mixing_bijections, NULL);
	make_block_invertible_matrix_pair(&final_encoding, NULL, 128);
	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand();
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);
	make_sbox_pair_8(input_sbox, input_sbox_inv);
	make_sbox_pair_128(middle_sbox, middle_sbox_inv);
	make_sbox_pair_8(output_sbox, output_sbox_inv);

	make_typeIB(typeIBs, tbox[NR - 1], final_encoding, input_sbox_inv,
			middle_sbox);
	make_typeIV_IB(typeIV_Is, middle_sbox_inv, output_sbox, output_sbox_inv);

	/*  prepare encoded input to the step */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			pre_enc_state[i][j] = rand();
			state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
					input_sbox[i][j]);
		}
	}
	shift_rows(state);
	do_typeIB(strips, state, typeIBs);
	do_typeIV_IB(state, strips, typeIV_Is);
	/*  decode output and present */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			state[i][j] = sub_bytes_hi_lo(state[i][j],
					output_sbox_inv[i][j]);
		}
	}
	dump_state("state: ", state);

	/*  equivalent transformation */
	slice_matrix_vertically(encoding_slices, 4 * 4, final_encoding);
	shift_rows(pre_enc_state);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row * 4 + col;
			uint8_t tbox_val = tbox[NR - 1][row][col][pre_enc_state[row][col]];
			mul_byte_by_matrix_128x8(strips[row][col],
					encoding_slices[index],
					tbox_val);
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			state2[row][col] = strips[15/4][15%4][index];
			for (k = 14; k != -1; --k) {
				state2[row][col] = strips[k/4][k%4][index] ^ state2[row][col];
			}
		}
	}
	dump_state("comparative state: ", state2);
	ASSERT(comp_states(state, state2));

	free_matrices(encoding_slices, 4*4);
	free_matrix(final_encoding);
	free_tbox_mixing_bijections(tbox_mixing_bijections);
}
void test_typeII()
{
	gf2matrix *mix_columns_mixing_bijection;
	tbox_mixing_bijections_t tbox_mixing_bijections;
	typeII_t typeIIs;
	tbox_t tbox;
	uint32_t expanded_key[(NR + 1) * 4];
	sboxes_8bit_t input_sbox[NR - 1], input_sbox_inv[NR - 1];
	sboxes_32bit_t output_sbox[NR - 1], output_sbox_inv[NR - 1];
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	_4bit_strip32_t strips;
	_4bit_strip32_t strips2;
	int round, row, col, i, j;

	make_block_invertible_matrix_pair(&mix_columns_mixing_bijection, NULL, 32);
	make_tbox_mixing_bijections(tbox_mixing_bijections, NULL);
	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand(); /*  identity key expansion */
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);
	make_rounds_sbox_pair_8(input_sbox, input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(output_sbox, output_sbox_inv, NR - 1);

	make_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
			input_sbox_inv, output_sbox);

	for (round = 0; round < NR - 1; ++round) {
		printf("round %d\n", round);
		/*  prepare encoded input to the step */
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				pre_enc_state[i][j] = rand();
				state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
						input_sbox[round][i][j]);
			}
		}
		shift_rows(state); /*  there's always shift rows before Type II */
		do_typeII(strips, state, typeIIs[round]);
		/*  decode output and present */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				sub_bytes_hi_lo32(strips[row][col],
						output_sbox_inv[round][row][col]);
			}
		}
		dump_4bit_strip32("strips: ", strips);

		/*  comparative section */
		shift_rows(pre_enc_state); /*  there's always shift rows before Type II */
		for (col = 0; col < 4; ++col) {
			uint8_t slice[4];
			for (row = 0; row < 4; ++row) {
				int new_row;
				uint8_t temp = tbox[round][row][col][pre_enc_state[row][col]];
				calc_mixing_slices(slice, temp, row);
				mul_array_by_matrix_32x32(slice,
						mix_columns_mixing_bijection,
						slice);
				for (new_row = 0; new_row < 4; ++new_row) {
					strips2[new_row][col][row] = slice[new_row];
				}
			}
		}
		dump_4bit_strip32("Comparative strips: ", strips2);
		ASSERT(comp_4bit_strips32(strips, strips2));
	}
	free_matrix(mix_columns_mixing_bijection);
	free_tbox_mixing_bijections(tbox_mixing_bijections);
}
void test_typeII_IV_combined()
{
	gf2matrix *mix_columns_mixing_bijection;
	tbox_mixing_bijections_t tbox_mixing_bijections;
	uint8_t typeIIs[NR - 1][4][4][256][4];
	uint8_t typeIV_IIs[NR - 1][4][4][3][2][128];
	tbox_t tbox;
	uint32_t expanded_key[(NR + 1) * 4];
	sboxes_8bit_t input_sbox[NR - 1], input_sbox_inv[NR - 1];
	sboxes_32bit_t middle_sbox[NR - 1], middle_sbox_inv[NR - 1];
	sboxes_8bit_t output_sbox[NR - 1], output_sbox_inv[NR - 1];
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	_4bit_strip32_t strips;
	_4bit_strip32_t strips2;
	int round, row, col, i, j, k;
	uint8_t state2[4][4];

	make_block_invertible_matrix_pair(&mix_columns_mixing_bijection, NULL, 32);
	/* &inv_mix_columns_mixing_bijection, 32); */
	make_tbox_mixing_bijections(tbox_mixing_bijections, NULL);
	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand(); /*  identity key expansion */
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);
	make_rounds_sbox_pair_8(input_sbox, input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(middle_sbox, middle_sbox_inv, NR - 1);
	make_rounds_sbox_pair_8(output_sbox, output_sbox_inv, NR - 1);

	make_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
			input_sbox_inv, middle_sbox);
	make_typeIV_II(typeIV_IIs, middle_sbox_inv, output_sbox, output_sbox_inv);

	for (round = 0; round < NR - 1; ++round) {
		printf("round %d\n", round);
		/*  prepare encoded input to the step */
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				pre_enc_state[i][j] = rand();
				state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
						input_sbox[round][i][j]);
			}
		}
		/* 		dump_state("pre_enc_state: ", pre_enc_state); */
		/* 		dump_state("state before: ", state); */
		shift_rows(state);
		do_typeII(strips, state, typeIIs[round]);
		do_typeIV_II(state, strips, typeIV_IIs[round]);
		/*  additionally decode the output and present the state */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = sub_bytes_hi_lo(state[row][col],
						output_sbox_inv[round][row][col]);
			}
		}
		dump_state("state: ", state);

		/*  comparative section */
		shift_rows(pre_enc_state);
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				uint8_t temp = tbox[round][row][col][pre_enc_state[row][col]];
				uint8_t slice[4];
				calc_mixing_slices(slice, temp, row);
				mul_array_by_matrix_32x32(slice,
						mix_columns_mixing_bijection,
						slice);
				for (k = 0; k < 4; ++k) {
					strips2[k][col][row] = slice[k];
				}
			}
		}
		/* dump_4bit_strip32("comp 4bitStrips: ", strips2); */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state2[row][col] = strips2[row][col][0] ^ strips2[row][col][1]
						^ strips2[row][col][2] ^ strips2[row][col][3];
			}
		}
		dump_state("Comparative state: ", state2);
		ASSERT(comp_states(state, state2));
	}
	free_matrix(mix_columns_mixing_bijection);
	free_tbox_mixing_bijections(tbox_mixing_bijections);
}
void test_typeIII()
{
	gf2matrix *inv_mix_columns_mixing_bijection;
	tbox_mixing_bijections_t inv_tbox_mixing_bijections;
	typeIII_t typeIIIs;
	sboxes_8bit_t input_sbox[NR - 1], input_sbox_inv[NR - 1];
	sboxes_32bit_t output_sbox[NR - 1], output_sbox_inv[NR - 1];
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	_4bit_strip32_t strips;
	_4bit_strip32_t strips2;
	int round, row, col, i, j, k;

	make_block_invertible_matrix_pair(&inv_mix_columns_mixing_bijection, NULL,
			32);
	make_tbox_mixing_bijections(inv_tbox_mixing_bijections, NULL);
	make_rounds_sbox_pair_8(input_sbox, input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(output_sbox, output_sbox_inv, NR - 1);
	make_typeIII(typeIIIs, inv_mix_columns_mixing_bijection,
			inv_tbox_mixing_bijections, input_sbox_inv, output_sbox);

	for (round = 0; round < NR - 1; ++round) {
		printf("round %d\n", round);
		/*  prepare encoded input to the step */
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				pre_enc_state[i][j] = rand();
				state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
						input_sbox[round][i][j]);
			}
		}
		do_typeIII(strips, state, typeIIIs[round]);
		/*  decode output and present */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				sub_bytes_hi_lo32(strips[row][col],
						output_sbox_inv[round][row][col]);
			}
		}
		dump_4bit_strip32("4bitstrips: ", strips);

		/*  comparative section */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				int new_row;
				uint8_t mc[4] = { 0, 0, 0, 0 };
				mc[row] = pre_enc_state[row][col];
				mul_array_by_matrix_32x32(mc,
						inv_mix_columns_mixing_bijection,
						mc);
				for (k = 0; k < 4; ++k) {
					mul_byte_by_matrix(&mc[k],
							inv_tbox_mixing_bijections[round][k][col],
							mc[k]);
				}
				for (new_row = 0; new_row < 4; ++new_row) {
					strips2[new_row][col][row] = mc[new_row];
				}
			}
		}
		dump_4bit_strip32("comparative 4bitstrips: ", strips2);
		ASSERT(comp_4bit_strips32(strips, strips2));
	}
	free_matrix(inv_mix_columns_mixing_bijection);
	free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
}

void test_typeIII_IV_combined()
{
	gf2matrix *inv_mix_columns_mixing_bijection;
	tbox_mixing_bijections_t inv_tbox_mixing_bijections;
	typeIII_t typeIIIs;
	typeIV_III_round_t typeIV_IIIs[NR - 1];
	sboxes_8bit_t input_sbox[NR - 1], input_sbox_inv[NR - 1];
	sboxes_32bit_t middle_sbox[NR - 1], middle_sbox_inv[NR - 1];
	sboxes_8bit_t output_sbox[NR - 1], output_sbox_inv[NR - 1];
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	_4bit_strip32_t strips;
	_4bit_strip32_t strips2;
	int round, row, col, i, j, k;
	uint8_t state2[4][4];

	make_block_invertible_matrix_pair(&inv_mix_columns_mixing_bijection, NULL,
			32);
	make_tbox_mixing_bijections(inv_tbox_mixing_bijections, NULL);
	make_rounds_sbox_pair_8(input_sbox, input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(middle_sbox, middle_sbox_inv, NR - 1);
	make_rounds_sbox_pair_8(output_sbox, output_sbox_inv, NR - 1);

	make_typeIII(typeIIIs, inv_mix_columns_mixing_bijection,
			inv_tbox_mixing_bijections, input_sbox_inv, middle_sbox);
	make_typeIV_III(typeIV_IIIs, middle_sbox_inv, output_sbox, output_sbox_inv);

	for (round = 0; round < NR - 1; ++round) {
		printf("round %d\n", round);
		/*  prepare encoded input to the step */
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				pre_enc_state[i][j] = rand();
				state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
						input_sbox[round][i][j]);
			}
		}
		/* 		dump_state("pre_enc_state: ", pre_enc_state); */
		/* 		dump_state("state before: ", state); */
		do_typeIII(strips, state, typeIIIs[round]);
		do_typeIV_III(state, strips, typeIV_IIIs[round]);
		/*  additionally decode the output and present the state */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = sub_bytes_hi_lo(state[row][col],
						output_sbox_inv[round][row][col]);
			}
		}
		dump_state("state: ", state);

		/*  comparative section */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				int new_row;
				uint8_t mc[4] = { 0, 0, 0, 0 };
				mc[row] = pre_enc_state[row][col];
				mul_array_by_matrix_32x32(mc,
						inv_mix_columns_mixing_bijection,
						mc);
				for (k = 0; k < 4; ++k) {
					mul_byte_by_matrix(&mc[k],
							inv_tbox_mixing_bijections[round][k][col],
							mc[k]);
				}
				for (new_row = 0; new_row < 4; ++new_row) {
					strips2[new_row][col][row] = mc[new_row];
				}
			}
		}
		/* dump_4bit_strip32("comp 4bitStrips: ", strips2); */
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state2[row][col] = strips2[row][col][0] ^ strips2[row][col][1]
						^ strips2[row][col][2] ^ strips2[row][col][3];
			}
		}
		dump_state("Comparative state: ", state2);
		ASSERT(comp_states(state, state2));
	}
	free_matrix(inv_mix_columns_mixing_bijection);
	free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
}
void test_typeIA_II_combined()
{
	gf2matrix *mix_columns_mixing_bijection;
	gf2matrix *initial_decoding;
	tbox_mixing_bijections_t tbox_mixing_bijections, inv_tbox_mixing_bijections;
	typeIA_t typeIAs;
	typeIV_IA_t typeIV_Is;
	typeII_t typeIIs;
	typeIV_II_round_t typeIV_IIs[NR - 1];
	tbox_t tbox;
	uint32_t expanded_key[(NR + 1) * 4];
	sboxes_8bit_t typeIA_input_sbox, typeIA_input_sbox_inv;
	sboxes_128bit_t typeIA_interim_sbox, typeIA_interim_sbox_inv;
	sboxes_8bit_t typeII_input_sbox[NR - 1], typeII_input_sbox_inv[NR - 1];
	sboxes_32bit_t typeII_interim_sbox[NR - 1], typeII_interim_sbox_inv[NR - 1];
	sboxes_8bit_t typeII_output_sbox[NR - 1], typeII_output_sbox_inv[NR - 1];
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	_4bit_strip32_t strips32;
	_4bit_strip128_t strips128;
	uint8_t state2[4][4];
	gf2matrix *decoding_slices[4 * 4];
	int round, row, col, i, j, k;

	make_block_invertible_matrix_pair(&mix_columns_mixing_bijection, NULL, 32);
	make_tbox_mixing_bijections(tbox_mixing_bijections,
			inv_tbox_mixing_bijections);
	make_block_invertible_matrix_pair(&initial_decoding, NULL, 128);
	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand(); /*  identity key expansion */
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);

	make_sbox_pair_8(typeIA_input_sbox, typeIA_input_sbox_inv);
	make_sbox_pair_128(typeIA_interim_sbox, typeIA_interim_sbox_inv);
	make_rounds_sbox_pair_8(typeII_input_sbox, typeII_input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(typeII_interim_sbox, typeII_interim_sbox_inv, NR
			- 1);
	make_rounds_sbox_pair_8(typeII_output_sbox, typeII_output_sbox_inv, NR - 1);

	make_typeIA(typeIAs, inv_tbox_mixing_bijections[0], initial_decoding,
			typeIA_input_sbox_inv, typeIA_interim_sbox);
	make_typeIV_IA(typeIV_Is, typeIA_interim_sbox_inv, typeII_input_sbox[0],
			typeII_input_sbox_inv[0]);

	make_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
			typeII_input_sbox_inv, typeII_interim_sbox);
	make_typeIV_II(typeIV_IIs, typeII_interim_sbox_inv, typeII_output_sbox,
			typeII_output_sbox_inv);

	round = 0;
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			pre_enc_state[i][j] = rand();
			state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
					typeIA_input_sbox[i][j]);
		}
	}
	/* 	dump_state("intial state: ", state); */
	/* 	dump_state("pre_enc_state: ", pre_enc_state); */

	do_typeIA(strips128, state, typeIAs);
	do_typeIV_IA(state, strips128, typeIV_Is);
	shift_rows(state);
	do_typeII(strips32, state, typeIIs[round]);
	do_typeIV_II(state, strips32, typeIV_IIs[round]);
	/*  additionally decode the output and present the state */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			state[row][col] = sub_bytes_hi_lo(state[row][col],
					typeII_output_sbox_inv[round][row][col]);
		}
	}
	dump_state("state: ", state);

	/*  comparative section */
	slice_matrix_vertically(decoding_slices, 4 * 4, initial_decoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			mul_byte_by_matrix_128x8(strips128[row][col], decoding_slices[row
					* 4 + col], pre_enc_state[row][col]);
			for (k = 0; k < 16; ++k) {
				mul_byte_by_matrix(&strips128[row][col][k],
						inv_tbox_mixing_bijections[0][k/4][k%4],
						strips128[row][col][k]);
			}
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			state2[row][col] = strips128[15/4][15%4][index];
			for (k = 14; k != -1; --k) {
				state2[row][col] = strips128[k/4][k%4][index] ^ state2[row][col];
			}
		}
	}

	shift_rows(state2);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			uint8_t temp = tbox[round][row][col][state2[row][col]];
			uint8_t slice[4];
			calc_mixing_slices(slice, temp, row);
			mul_array_by_matrix_32x32(slice,
					mix_columns_mixing_bijection,
					slice);
			for (k = 0; k < 4; ++k) {
				strips32[k][col][row] = slice[k];
			}
		}
	}
	/* dump_4bit_strip32("comp 4bitStrips: ", strips2); */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			state2[row][col] = strips32[row][col][0] ^ strips32[row][col][1]
					^ strips32[row][col][2] ^ strips32[row][col][3];
		}
	}
	dump_state("Comparative state: ", state2);
	ASSERT(comp_states(state, state2));
	free_matrices(decoding_slices, 4 * 4);
	free_matrix(initial_decoding);
	free_tbox_mixing_bijections(tbox_mixing_bijections);
	free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
	free_matrix(mix_columns_mixing_bijection);
}
void test_typeIII_IB_combined()
{
	gf2matrix *inv_mix_columns_mixing_bijection;
	tbox_mixing_bijections_t tbox_mixing_bijections, inv_tbox_mixing_bijections;
	gf2matrix *final_encoding;
	uint32_t expanded_key[(NR + 1) * 4];
	tbox_t tbox;
	typeIB_t typeIBs;
	typeIV_IB_t typeIV_IBs;
	typeII_t typeIIIs;
	typeIV_III_round_t typeIV_IIIs[NR - 1];
	sboxes_8bit_t typeIII_input_sbox[NR - 1], typeIII_input_sbox_inv[NR - 1];
	sboxes_32bit_t typeIII_middle_sbox[NR - 1], typeIII_middle_sbox_inv[NR - 1];
	sboxes_8bit_t typeIII_output_sbox[NR - 1], typeIII_output_sbox_inv[NR - 1];
	sboxes_128bit_t typeIB_middle_sbox, typeIB_middle_sbox_inv;
	sboxes_8bit_t typeIB_output_sbox, typeIB_output_sbox_inv;
	_4bit_strip128_t strips128;
	_4bit_strip32_t strips32;
	uint8_t state[4][4];
	uint8_t pre_enc_state[4][4];
	int round, row, col, i, j, k;

	make_block_invertible_matrix_pair(&inv_mix_columns_mixing_bijection, NULL,
			32);
	make_tbox_mixing_bijections(tbox_mixing_bijections,
			inv_tbox_mixing_bijections);
	make_block_invertible_matrix_pair(&final_encoding, NULL, 128);
	for (i = 0; i < (NR + 1) * 4; ++i)
		expanded_key[i] = (uint8_t) rand(); /*  identity key expansion */
	make_tbox(tbox, SBox, expanded_key, tbox_mixing_bijections);

	make_rounds_sbox_pair_8(typeIII_input_sbox, typeIII_input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(typeIII_middle_sbox, typeIII_middle_sbox_inv, NR
			- 1);
	make_rounds_sbox_pair_8(typeIII_output_sbox, typeIII_output_sbox_inv, NR
			- 1);
	make_sbox_pair_128(typeIB_middle_sbox, typeIB_middle_sbox_inv);
	make_sbox_pair_8(typeIB_output_sbox, typeIB_output_sbox_inv);

	make_typeIII(typeIIIs, inv_mix_columns_mixing_bijection,
			inv_tbox_mixing_bijections, typeIII_input_sbox_inv,
			typeIII_middle_sbox);
	make_typeIV_III(typeIV_IIIs, typeIII_middle_sbox_inv, typeIII_output_sbox,
			typeIII_output_sbox_inv);
	make_typeIB(typeIBs, tbox[NR - 1], final_encoding,
			typeIII_output_sbox_inv[NR - 2], typeIB_middle_sbox);
	make_typeIV_IB(typeIV_IBs, typeIB_middle_sbox_inv, typeIB_output_sbox,
			typeIB_output_sbox_inv);

	round = NR - 2;

	/*  prepare encoded input to the step */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			pre_enc_state[i][j] = rand();
			state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
					typeIII_input_sbox[round][i][j]);
		}
	}
	do_typeIII(strips32, state, typeIIIs[round]);
	do_typeIV_III(state, strips32, typeIV_IIIs[round]);
	shift_rows(state);
	do_typeIB(strips128, state, typeIBs);
	do_typeIV_IB(state, strips128, typeIV_IBs);
	/*  decode output and present */
	for (i = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j) {
			state[i][j] = sub_bytes_hi_lo(state[i][j],
					typeIB_output_sbox_inv[i][j]);
		}
	}
	dump_state("state: ", state);

	/*  comparative section */
	uint8_t state2[4][4];
	gf2matrix *encoding_slices[4 * 4];
	slice_matrix_vertically(encoding_slices, 4 * 4, final_encoding);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int new_row;
			uint8_t mc[4] = { 0, 0, 0, 0 };
			mc[row] = pre_enc_state[row][col];
			mul_array_by_matrix_32x32(mc, inv_mix_columns_mixing_bijection, mc);
			for (k = 0; k < 4; ++k) {
				mul_byte_by_matrix(&mc[k],
						inv_tbox_mixing_bijections[round][k][col],
						mc[k]);
			}
			for (new_row = 0; new_row < 4; ++new_row) {
				strips32[new_row][col][row] = mc[new_row];
			}
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			state2[row][col] = strips32[row][col][0] ^ strips32[row][col][1]
					^ strips32[row][col][2] ^ strips32[row][col][3];
		}
	}
	shift_rows(state2);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row * 4 + col;
			uint8_t tbox_val = tbox[NR - 1][row][col][state2[row][col]];
			mul_byte_by_matrix_128x8(strips128[row][col],
					encoding_slices[index],
					tbox_val);
		}
	}
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			state2[row][col] = strips128[15/4][15%4][index];
			for (k = 14; k != -1; --k) {
				state2[row][col] = strips128[k/4][k%4][index] ^ state2[row][col];
			}
		}
	}

	dump_state("comparative state: ", state2);
	ASSERT(comp_states(state, state2));

	free_matrices(encoding_slices, 4*4);
	free_matrix(inv_mix_columns_mixing_bijection);
	free_tbox_mixing_bijections(tbox_mixing_bijections);
	free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
	free_matrix(final_encoding);
}
void test_encryption()
{
	gf2matrix *mix_columns_mixing_bijection, *inv_mix_columns_mixing_bijection;
	tbox_mixing_bijections_t tbox_mixing_bijections, inv_tbox_mixing_bijections;
	gf2matrix *initial_encoding, *initial_decoding;
	gf2matrix *final_encoding, *final_decoding;
	tbox_t tbox;
	typeIA_t typeIAs;
	typeII_t typeIIs;
	typeIII_t typeIIIs;
	typeIB_t typeIBs;
	typeIV_IA_t typeIV_IAs;
	typeIV_IB_t typeIV_IBs;
	typeIV_II_round_t typeIV_IIs[NR - 1];
	typeIV_III_round_t typeIV_IIIs[NR - 1];
	uint32_t key_schedule[4 * (NR + 1)];
	uint8_t key[KEY_SIZE];

	sboxes_8bit_t typeIA_input_sbox, typeIA_input_sbox_inv;
	sboxes_128bit_t typeIA_interim_sbox, typeIA_interim_sbox_inv;
	sboxes_8bit_t typeII_input_sbox[NR], typeII_input_sbox_inv[NR];
	sboxes_32bit_t typeII_interim_sbox[NR - 1], typeII_interim_sbox_inv[NR - 1];
	sboxes_8bit_t typeII_output_sbox[NR - 1], typeII_output_sbox_inv[NR - 1];
	sboxes_32bit_t typeIII_interim_sbox[NR - 1], typeIII_interim_sbox_inv[NR
			- 1];
	sboxes_128bit_t typeIB_interim_sbox, typeIB_interim_sbox_inv;
	sboxes_8bit_t typeIB_output_sbox, typeIB_output_sbox_inv;
	uint8_t state[4][4];
	uint8_t in[16];
	uint8_t out[16], out2[16];
	_4bit_strip32_t strips32;
	_4bit_strip128_t strips128;
	int round, row, col, i;

	int tries = 3;
	for (; tries != 0; --tries) {
		randomize_key(key);
		make_block_invertible_matrix_pair(&mix_columns_mixing_bijection,
				&inv_mix_columns_mixing_bijection, 32);
		make_tbox_mixing_bijections(tbox_mixing_bijections,
				inv_tbox_mixing_bijections);

		make_block_invertible_matrix_pair(&initial_encoding, &initial_decoding, 128);
		make_block_invertible_matrix_pair(&final_encoding, &final_decoding, 128);
		expand_key(key, SBox, key_schedule, 4);
		make_tbox(tbox, SBox, key_schedule, tbox_mixing_bijections);

		make_sbox_pair_8(typeIA_input_sbox, typeIA_input_sbox_inv);
		make_sbox_pair_128(typeIA_interim_sbox, typeIA_interim_sbox_inv);
		make_rounds_sbox_pair_8(typeII_input_sbox, typeII_input_sbox_inv, NR);
		make_rounds_sbox_pair_32(typeII_interim_sbox, typeII_interim_sbox_inv,
				NR - 1);
		make_rounds_sbox_pair_8(typeII_output_sbox, typeII_output_sbox_inv, NR - 1);
		make_rounds_sbox_pair_32(typeIII_interim_sbox, typeIII_interim_sbox_inv,
				NR - 1);
		make_sbox_pair_128(typeIB_interim_sbox, typeIB_interim_sbox_inv);
		make_sbox_pair_8(typeIB_output_sbox, typeIB_output_sbox_inv);

		make_typeIA(typeIAs, inv_tbox_mixing_bijections[0], initial_decoding,
				typeIA_input_sbox_inv, typeIA_interim_sbox);
		make_typeIV_IA(typeIV_IAs, typeIA_interim_sbox_inv, typeII_input_sbox[0],
				typeII_input_sbox_inv[0]);
		make_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
				typeII_input_sbox_inv, typeII_interim_sbox);
		make_typeIV_II(typeIV_IIs, typeII_interim_sbox_inv, typeII_output_sbox,
				typeII_output_sbox_inv);
		make_typeIII(typeIIIs, inv_mix_columns_mixing_bijection,
				&inv_tbox_mixing_bijections[1], typeII_output_sbox_inv,
				typeIII_interim_sbox);
		make_typeIV_III(typeIV_IIIs, typeIII_interim_sbox_inv,
				&typeII_input_sbox[1], &typeII_input_sbox_inv[1]);
		make_typeIB(typeIBs, tbox[NR - 1], final_encoding,
				typeII_input_sbox_inv[NR - 1], typeIB_interim_sbox);
		make_typeIV_IB(typeIV_IBs, typeIB_interim_sbox_inv, typeIB_output_sbox,
				typeIB_output_sbox_inv);

		for (i = 0; i < 16; ++i) {
			in[i] = rand();
		}
		dump_hex("input: ", in, 16);
		printf("White-box AES Cipher:\n");
		do_input(state, in, initial_encoding, typeIA_input_sbox);

		dump_state("State before ", state);
		do_typeIA(strips128, state, typeIAs);
		do_typeIV_IA(state, strips128, typeIV_IAs);
		for (round = 0; round < NR - 1; ++round) {
			printf("round %d: ", round);
			dump_state("", state);
			shift_rows(state);
			do_typeII(strips32, state, typeIIs[round]);
			do_typeIV_II(state, strips32, typeIV_IIs[round]);
			do_typeIII(strips32, state, typeIIIs[round]);
			do_typeIV_III(state, strips32, typeIV_IIIs[round]);
		}
		shift_rows(state);
		dump_state("round 9-10: ", state);
		do_typeIB(strips128, state, typeIBs);
		do_typeIV_IB(state, strips128, typeIV_IBs);

		do_output(out, state, final_decoding, typeIB_output_sbox_inv);

		printf("Original AES Cipher on same input using same key:\n");
		cipher(in, out2, SBox, key_schedule);
		dump_hex("WB Output ", out, 16);
		dump_hex("AES Output ", out2, 16);
		ASSERT(memcmp(out, out2, 16) == 0);
		free_matrix(mix_columns_mixing_bijection);
		free_matrix(inv_mix_columns_mixing_bijection);
		free_tbox_mixing_bijections(tbox_mixing_bijections);
		free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
		free_matrix(final_decoding);
		free_matrix(final_encoding);
		free_matrix(initial_decoding);
		free_matrix(initial_encoding);
	}
}
void test_encryptdecrypt()
{
	gf2matrix *mix_columns_mixing_bijection, *inv_mix_columns_mixing_bijection;
	tbox_mixing_bijections_t tbox_mixing_bijections, inv_tbox_mixing_bijections;
	gf2matrix *initial_encoding, *initial_decoding;
	gf2matrix *final_encoding, *final_decoding;
	tbox_t tbox;
	typeIA_t typeIAs;
	typeII_t typeIIs;
	typeIII_t typeIIIs;
	typeIB_t typeIBs;
	typeIV_IA_t typeIV_IAs;
	typeIV_IB_t typeIV_IBs;
	typeIV_II_round_t typeIV_IIs[NR - 1];
	typeIV_III_round_t typeIV_IIIs[NR - 1];
	uint32_t key_schedule[4 * (NR + 1)];
	uint8_t key[KEY_SIZE];

	sboxes_8bit_t typeIA_input_sbox, typeIA_input_sbox_inv;
	sboxes_128bit_t typeIA_interim_sbox, typeIA_interim_sbox_inv;
	sboxes_8bit_t typeII_input_sbox[NR], typeII_input_sbox_inv[NR];
	sboxes_32bit_t typeII_interim_sbox[NR - 1], typeII_interim_sbox_inv[NR - 1];
	sboxes_8bit_t typeII_output_sbox[NR - 1], typeII_output_sbox_inv[NR - 1];
	sboxes_32bit_t typeIII_interim_sbox[NR - 1], typeIII_interim_sbox_inv[NR
			- 1];
	sboxes_128bit_t typeIB_interim_sbox, typeIB_interim_sbox_inv;
	sboxes_8bit_t typeIB_output_sbox, typeIB_output_sbox_inv;
	uint8_t state[4][4];
	uint8_t in[16];
	uint8_t out[16], out2[16];
	_4bit_strip32_t strips32;
	_4bit_strip128_t strips128;
	int round, row, col, i;
	uint32_t mixed_key_schedule[4 * (NR + 1)];
	uint8_t dec_out[16];

	int tries = 3;
	for (; tries != 0; --tries) {
		randomize_key(key);
		make_block_invertible_matrix_pair(&mix_columns_mixing_bijection,
				&inv_mix_columns_mixing_bijection, 32);
		make_tbox_mixing_bijections(tbox_mixing_bijections,
				inv_tbox_mixing_bijections);

		make_block_invertible_matrix_pair(&initial_encoding, &initial_decoding, 128);
		make_block_invertible_matrix_pair(&final_encoding, &final_decoding, 128);

		expand_key(key, SBox, key_schedule, 4);
		make_tbox(tbox, SBox, key_schedule, tbox_mixing_bijections);

		make_sbox_pair_8(typeIA_input_sbox, typeIA_input_sbox_inv);
		make_sbox_pair_128(typeIA_interim_sbox, typeIA_interim_sbox_inv);
		make_rounds_sbox_pair_8(typeII_input_sbox, typeII_input_sbox_inv, NR);
		make_rounds_sbox_pair_32(typeII_interim_sbox, typeII_interim_sbox_inv,
				NR - 1);
		make_rounds_sbox_pair_8(typeII_output_sbox, typeII_output_sbox_inv, NR - 1);
		make_rounds_sbox_pair_32(typeIII_interim_sbox, typeIII_interim_sbox_inv,
				NR - 1);
		make_sbox_pair_128(typeIB_interim_sbox, typeIB_interim_sbox_inv);
		make_sbox_pair_8(typeIB_output_sbox, typeIB_output_sbox_inv);

		make_typeIA(typeIAs, inv_tbox_mixing_bijections[0], initial_decoding,
				typeIA_input_sbox_inv, typeIA_interim_sbox);
		make_typeIV_IA(typeIV_IAs, typeIA_interim_sbox_inv, typeII_input_sbox[0],
				typeII_input_sbox_inv[0]);
		make_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
				typeII_input_sbox_inv, typeII_interim_sbox);
		make_typeIV_II(typeIV_IIs, typeII_interim_sbox_inv, typeII_output_sbox,
				typeII_output_sbox_inv);
		make_typeIII(typeIIIs, inv_mix_columns_mixing_bijection,
				&inv_tbox_mixing_bijections[1], typeII_output_sbox_inv,
				typeIII_interim_sbox);
		make_typeIV_III(typeIV_IIIs, typeIII_interim_sbox_inv,
				&typeII_input_sbox[1], &typeII_input_sbox_inv[1]);
		make_typeIB(typeIBs, tbox[NR - 1], final_encoding,
				typeII_input_sbox_inv[NR - 1], typeIB_interim_sbox);
		make_typeIV_IB(typeIV_IBs, typeIB_interim_sbox_inv, typeIB_output_sbox,
				typeIB_output_sbox_inv);

		for (i = 0; i < 16; ++i) {
			in[i] = rand();
		}
		dump_hex("input: ", in, 16);
		printf("White-box AES Cipher:\n");
		do_input(state, in, initial_encoding, typeIA_input_sbox);

		dump_state("State before ", state);
		do_typeIA(strips128, state, typeIAs);
		do_typeIV_IA(state, strips128, typeIV_IAs);
		for (round = 0; round < NR - 1; ++round) {
			printf("round %d: ", round);
			dump_state("", state);
			shift_rows(state);
			do_typeII(strips32, state, typeIIs[round]);
			do_typeIV_II(state, strips32, typeIV_IIs[round]);
			do_typeIII(strips32, state, typeIIIs[round]);
			do_typeIV_III(state, strips32, typeIV_IIIs[round]);
		}
		shift_rows(state);
		dump_state("round 9-10: ", state);
		do_typeIB(strips128, state, typeIBs);
		do_typeIV_IB(state, strips128, typeIV_IBs);

		do_output(out, state, final_decoding, typeIB_output_sbox_inv);
		dump_hex("WB Output ", out, 16);

		printf("Original AES Cipher on same input using same key:\n");
		cipher(in, out2, SBox, key_schedule);
		dump_hex("AES Output ", out2, 16);
		ASSERT(memcmp(out, out2, 16) == 0);
		{
			uint8_t decipher_out[4*4];
			printf("Comparing AES decipher with the input\n");
			decipher(out, decipher_out, ISBox, key_schedule);
			dump_hex("AES decipher Output ", decipher_out, 16);
			ASSERT(memcmp(in, decipher_out, 16) == 0);
		}

		// start decrypting the same thing
		// at the end we should get the input
		make_block_invertible_matrix_pair(&mix_columns_mixing_bijection,
				&inv_mix_columns_mixing_bijection, 32);
		make_tbox_mixing_bijections(tbox_mixing_bijections,
				inv_tbox_mixing_bijections);

		make_block_invertible_matrix_pair(&initial_encoding, &initial_decoding, 128);
		make_block_invertible_matrix_pair(&final_encoding, &final_decoding, 128);

		//expand_key(key, SBox, mixed_key_schedule, 4);
		mix_expanded_key(key_schedule);
		make_inv_tbox(tbox, ISBox, key_schedule, tbox_mixing_bijections);
		make_sbox_pair_8(typeIA_input_sbox, typeIA_input_sbox_inv);
		make_sbox_pair_128(typeIA_interim_sbox, typeIA_interim_sbox_inv);
		make_rounds_sbox_pair_8(typeII_input_sbox, typeII_input_sbox_inv, NR);
		make_rounds_sbox_pair_32(typeII_interim_sbox, typeII_interim_sbox_inv,
				NR - 1);
		make_rounds_sbox_pair_8(typeII_output_sbox, typeII_output_sbox_inv, NR - 1);
		make_rounds_sbox_pair_32(typeIII_interim_sbox, typeIII_interim_sbox_inv,
				NR - 1);
		make_sbox_pair_128(typeIB_interim_sbox, typeIB_interim_sbox_inv);
		make_sbox_pair_8(typeIB_output_sbox, typeIB_output_sbox_inv);


		make_typeIA(typeIAs, inv_tbox_mixing_bijections[NR - 1], initial_decoding,
				typeIA_input_sbox_inv, typeIA_interim_sbox);
		make_typeIV_IA(typeIV_IAs, typeIA_interim_sbox_inv,
				typeII_input_sbox[NR - 1], typeII_input_sbox_inv[NR-1]);
		make_inv_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
				&typeII_input_sbox_inv[1], typeII_interim_sbox);
		make_typeIV_II(typeIV_IIs, typeII_interim_sbox_inv, typeII_output_sbox,
				typeII_output_sbox_inv);
		make_typeIII(typeIIIs, inv_mix_columns_mixing_bijection,
				inv_tbox_mixing_bijections, typeII_output_sbox_inv,
				typeIII_interim_sbox);
		make_typeIV_III(typeIV_IIIs, typeIII_interim_sbox_inv,
				typeII_input_sbox, typeII_input_sbox_inv);
		make_inv_typeIB(typeIBs, tbox[0], final_encoding,
				typeII_input_sbox_inv[0], typeIB_interim_sbox);
		make_typeIV_IB(typeIV_IBs, typeIB_interim_sbox_inv, typeIB_output_sbox,
				typeIB_output_sbox_inv);

		// the input to this stage is the output of the encryption stage
		do_input(state, out, initial_encoding, typeIA_input_sbox);

		dump_state("State before ", state);
		do_typeIA(strips128, state, typeIAs);
		do_typeIV_IA(state, strips128, typeIV_IAs);
		for (round = NR - 2; round != -1; --round) {
			printf("round %d: ", round + 2);
			dump_state("", state);
			inv_shift_rows(state);
			do_typeII(strips32, state, typeIIs[round]);
			do_typeIV_II(state, strips32, typeIV_IIs[round]);
			do_typeIII(strips32, state, typeIIIs[round]);
			do_typeIV_III(state, strips32, typeIV_IIIs[round]);
		}
		inv_shift_rows(state);
		dump_state("rounds 1 and 0: ", state);
		do_typeIB(strips128, state, typeIBs);
		do_typeIV_IB(state, strips128, typeIV_IBs);

		do_output(dec_out, state, final_decoding, typeIB_output_sbox_inv);

		// compare the in and out
		dump_hex("WB-AES INV Output ", dec_out, 16);
		printf("Original AES Equivalent Inverse Cipher on same input using same key:\n");
		eqv_decipher(out, out2, ISBox, key_schedule);
		dump_hex("AES Output ", out2, 16);
		ASSERT(memcmp(out2, dec_out, 16) == 0);

		ASSERT(memcmp(dec_out, in, 16) == 0);
		// cleanup
		free_matrix(mix_columns_mixing_bijection);
		free_matrix(inv_mix_columns_mixing_bijection);
		free_tbox_mixing_bijections(tbox_mixing_bijections);
		free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
		free_matrix(final_decoding);
		free_matrix(final_encoding);
		free_matrix(initial_decoding);
		free_matrix(initial_encoding);
	}
}
int main()
{
	srandom(1234);
	TEST(test_xor_hi_lo);
	TEST(test_mix_columns_mixing_bijection);
	TEST(test_tbox);
	TEST(test_typeIA);
	TEST(test_typeIB);
	TEST(test_typeIV128);
	TEST(test_typeIV32);
	TEST(test_typeIA_IV_combined);
	TEST(test_typeIB_IV_combined);
	TEST(test_typeII);
	TEST(test_typeIII);
	TEST(test_typeII_IV_combined);
	TEST(test_typeIII_IV_combined);
	TEST(test_typeIA_II_combined);
	TEST(test_typeIII_IB_combined);
	TEST(test_encryption);
	TEST(test_encryptdecrypt);
	return 0;
}
