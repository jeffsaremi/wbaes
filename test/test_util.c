#include <memory.h>
#include "test_util.h"
#include "wb_aes.h"
#include "wb_aes_gen.h"

int comp_states(uint8_t state1[4][4], uint8_t state2[4][4])
{
	return memcmp(&state1[0][0], &state2[0][0], 16) == 0;
}
int comp_4bit_strips32(uint8_t _4bitStrip1[4][4][4],
		uint8_t _4bitStrip2[4][4][4])
{
	return memcmp(&_4bitStrip1[0][0][0], &_4bitStrip2[0][0][0], 16*4) == 0;
}
int comp_4bit_strips128(uint8_t _4bitStrip1[4][4][16],
		uint8_t _4bitStrip2[4][4][16])
{
	return memcmp(&_4bitStrip1[0][0][0], &_4bitStrip2[0][0][0], 16*16) == 0;
}
void copy_state_col(uint8_t col_bytes[4], uint8_t state[4][4], int col)
{
	int i;
	for (i = 0; i < 4; ++i)
		col_bytes[i] = state[i][col];
}
void copy_col2state(uint8_t state[4][4], uint8_t col_bytes[4], int col)
{
	int i;
	for (i = 0; i < 4; ++i)
		state[i][col] = col_bytes[i];
}
void randomize_state(uint8_t state[4][4])
{
	int i, j;
	for (i = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			state[i][j] = (uint8_t) rand();
}
void inv_mix_columns(uint8_t state[4][4])
{
	int col;
	for (col = 0; col < 4; ++col) {
		uint8_t slice0[4], slice1[4], slice2[4], slice3[4];
		calc_inv_mixing_slices(slice0, state[0][col], 0);
		calc_inv_mixing_slices(slice1, state[1][col], 1);
		calc_inv_mixing_slices(slice2, state[2][col], 2);
		calc_inv_mixing_slices(slice3, state[3][col], 3);

		state[0][col] = slice0[0] ^ slice1[0] ^ slice2[0] ^ slice3[0];
		state[1][col] = slice0[1] ^ slice1[1] ^ slice2[1] ^ slice3[1];
		state[2][col] = slice0[2] ^ slice1[2] ^ slice2[2] ^ slice3[2];
		state[3][col] = slice0[3] ^ slice1[3] ^ slice2[3] ^ slice3[3];
	}
}
void make_identity_tbox(tbox_t tbox)
{
	int round, i, j, x;
	for (round = 0; round < NR; ++round)
		for (i = 0; i < 4; ++i)
			for (j = 0; j < 4; ++j)
				for (x = 0; x < 256; ++x)
					tbox[round][i][j][x] = (uint8_t) x;
}
void make_identity_4bit_sbox(uint8_t sbox[16])
{
	int i;
	for (i = 0; i < 16; ++i)
		sbox[i] = i;
}
void make_identity_sboxes8(sboxes_8bit_t sboxes)
{
	int i, j, k;
	for (i = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			for (k = 0; k < 2; ++k)
				make_identity_4bit_sbox(sboxes[i][j][k]);
}
void make_identity_sboxes32(sboxes_32bit_t sboxes)
{
	int i, j, k;
	for (i = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			for (k = 0; k < 8; ++k)
				make_identity_4bit_sbox(sboxes[i][j][k]);
}
void make_identity_sboxes128(sboxes_128bit_t sboxes)
{
	int i, j, k;
	for (i = 0; i < 4; ++i)
		for (j = 0; j < 4; ++j)
			for (k = 0; k < 32; ++k)
				make_identity_4bit_sbox(sboxes[i][j][k]);
}
void make_identity_tbox_mixing_bijections(tbox_mixing_bijections_t tbox_mbs)
{
	int round, row, col;
	for (round = 0; round < NR; ++round) {
		for (row = 0; row < 4; ++row) {
			for (col = 0; col < 4; ++col) {
				tbox_mbs[round][row][col] = make_identity_matrix(8);
			}
		}
	}
}
void randomize_key(uint8_t *key)
{
	int i;
	for(i = 0; i < KEY_SIZE; ++i)
		key[i] = rand();
}
