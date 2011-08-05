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

void test_inv_tbox()
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
	make_inv_tbox(tbox, ISBox, expanded_key, tbox_mixing_bijections);

	for (round = NR - 1; round != 1; --round) {
		printf("round: %d \n", round);
		randomize_state(state);
		memcpy(&state2[0][0], &state[0][0], 16); /*  for validation */

		dump_state("State before ", state);
		inv_shift_rows(state);
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = tbox[round][row][col][(state[row][col])];
			}
		}
		dump_state("State after TBox ", state);
		/*  validation */
		/*  The above should be equal to: */
		/*  0. mix */
		/* 	1. addroundkey */
		/*  2. subbytes */
		/* 	3. shiftrows */

		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				mul_byte_by_matrix(&state2[row][col],
						tbox_mixing_bijections[round][row][col],
						state2[row][col]);
			}
		}
		add_round_key(state2, expanded_key, round+1);
		inv_sub_bytes(state2, ISBox);
		inv_shift_rows(state2);

		dump_state("Validation State ", state2);
		ASSERT(comp_states(state, state2));
	}
	/*  validation for the last round is different */
	round = 0;
	{
		printf("rounds 9 and 10\n");
		randomize_state(state);
		memcpy(&state2[0][0], &state[0][0], 16); /*  for validation */

		dump_state("State before ", state);
		inv_shift_rows(state);
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				state[row][col] = tbox[round][row][col][(state[row][col])];
			}
		}
		dump_state("State after TBox ", state);
		/*  validation */
		/*  The above should be equal to: */
		/*  0. mix */
		/* 	1. addroundkey */
		/*  2. subbytes */
		/* 	3. shiftrows */
		/* 	4. addroundkey */

		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				mul_byte_by_matrix(&state2[row][col],
						tbox_mixing_bijections[round][row][col],
						state2[row][col]);
			}
		}
		add_round_key(state2, expanded_key, 1);
		inv_sub_bytes(state2, ISBox);
		inv_shift_rows(state2);
		add_round_key(state2, expanded_key, 0); /*  the last key */

		dump_state("Validation State ", state2);
		ASSERT(comp_states(state, state2));
	}
	free_tbox_mixing_bijections(tbox_mixing_bijections);
	free_tbox_mixing_bijections(inv_tbox_mixing_bijections);
}
void test_inv_typeII()
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
	make_tbox(tbox, ISBox, expanded_key, tbox_mixing_bijections);
	make_rounds_sbox_pair_8(input_sbox, input_sbox_inv, NR - 1);
	make_rounds_sbox_pair_32(output_sbox, output_sbox_inv, NR - 1);

	make_inv_typeII(typeIIs, tbox, mix_columns_mixing_bijection,
			input_sbox_inv, output_sbox);

	for (round = NR - 2; round != -1; --round) {
		printf("round %d\n", round);
		/*  prepare encoded input to the step */
		for (i = 0; i < 4; ++i) {
			for (j = 0; j < 4; ++j) {
				pre_enc_state[i][j] = rand();
				state[i][j] = sub_bytes_hi_lo(pre_enc_state[i][j],
						input_sbox[round][i][j]);
			}
		}
		inv_shift_rows(state); /*  there's always shift rows before Type II */
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
		inv_shift_rows(pre_enc_state); /*  there's always shift rows before Type II */
		for (col = 0; col < 4; ++col) {
			uint8_t slice[4];
			for (row = 0; row < 4; ++row) {
				int new_row;
				uint8_t temp = tbox[round + 1][row][col][pre_enc_state[row][col]];
				calc_inv_mixing_slices(slice, temp, row);
				mul_array_by_matrix_32x32(slice, mix_columns_mixing_bijection,
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

void test_decryption()
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
	uint32_t mixed_key_schedule[4 * (NR + 1)];
	uint8_t key[KEY_SIZE] = {
			0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6,
			0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
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
/* 		mix_columns_mixing_bijection = make_identity_matrix(32); */
/* 		inv_mix_columns_mixing_bijection = make_identity_matrix(32); */
		make_tbox_mixing_bijections(tbox_mixing_bijections,
				inv_tbox_mixing_bijections);
/* 		make_identity_tbox_mixing_bijections(tbox_mixing_bijections); */
/* 		make_identity_tbox_mixing_bijections(inv_tbox_mixing_bijections); */

		make_block_invertible_matrix_pair(&initial_encoding, &initial_decoding, 128);
/* 		initial_encoding = make_identity_matrix(128); */
/* 		initial_decoding = make_identity_matrix(128); */
		make_block_invertible_matrix_pair(&final_encoding, &final_decoding, 128);
/* 		final_encoding = make_identity_matrix(128); */
/* 		final_decoding = make_identity_matrix(128); */

		expand_key(key, SBox, mixed_key_schedule, 4);
		mix_expanded_key(mixed_key_schedule);
		make_inv_tbox(tbox, ISBox, mixed_key_schedule, tbox_mixing_bijections);

		make_sbox_pair_8(typeIA_input_sbox, typeIA_input_sbox_inv);
/* 		make_identity_sboxes8(typeIA_input_sbox); */
/* 		make_identity_sboxes8(typeIA_input_sbox_inv); */
		make_sbox_pair_128(typeIA_interim_sbox, typeIA_interim_sbox_inv);
/* 		make_identity_sboxes128(typeIA_interim_sbox); */
/* 		make_identity_sboxes128(typeIA_interim_sbox_inv); */
		make_rounds_sbox_pair_8(typeII_input_sbox, typeII_input_sbox_inv, NR);
/* 		for(i =0; i < NR; ++i) { */
/* 			make_identity_sboxes8(typeII_input_sbox[i]); */
/* 			make_identity_sboxes8(typeII_input_sbox_inv[i]); */
/* 		} */
		make_rounds_sbox_pair_32(typeII_interim_sbox, typeII_interim_sbox_inv,
				NR - 1);
/* 		for(i =0; i < NR-1; ++i) { */
/* 			make_identity_sboxes32(typeII_interim_sbox[i]); */
/* 			make_identity_sboxes32(typeII_interim_sbox_inv[i]); */
/* 		} */
		make_rounds_sbox_pair_8(typeII_output_sbox, typeII_output_sbox_inv, NR - 1);
/* 		for(i =0; i < NR-1; ++i) { */
/* 			make_identity_sboxes8(typeII_output_sbox[i]); */
/* 			make_identity_sboxes8(typeII_output_sbox_inv[i]); */
/* 		} */

		make_rounds_sbox_pair_32(typeIII_interim_sbox, typeIII_interim_sbox_inv,
				NR - 1);
/* 		for(i =0; i < NR-1; ++i) { */
/* 			make_identity_sboxes32(typeIII_interim_sbox[i]); */
/* 			make_identity_sboxes32(typeIII_interim_sbox_inv[i]); */
/* 		} */
		make_sbox_pair_128(typeIB_interim_sbox, typeIB_interim_sbox_inv);
/* 		make_identity_sboxes128(typeIB_interim_sbox); */
/* 		make_identity_sboxes128(typeIB_interim_sbox_inv); */
		make_sbox_pair_8(typeIB_output_sbox, typeIB_output_sbox_inv);
/* 		make_identity_sboxes8(typeIB_output_sbox); */
/* 		make_identity_sboxes8(typeIB_output_sbox_inv); */


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

		for (i = 0; i < 16; ++i) {
			in[i] = rand();
		}
		dump_hex("input: ", in, 16);
		printf("White-box inverse cipher:\n");
		do_input(state, in, initial_encoding, typeIA_input_sbox);

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

		do_output(out, state, final_decoding, typeIB_output_sbox_inv);

		printf("Original AES Equivalent Inverse Cipher on same input using same key:\n");
		eqv_decipher(in, out2, ISBox, mixed_key_schedule);
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

int main()
{
	srandom(1234);
	TEST(test_inv_tbox);
	TEST(test_inv_typeII);
	TEST(test_decryption);
	return 0;
}
