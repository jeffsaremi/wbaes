#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include "wb_aes.h"
#include "util.h"

uint8_t sub_bytes_hi_lo(uint8_t x, uint8_t sbox[][16])
{
	return sbox[1][(x & 0x0F)] | (sbox[0][(x >> 4)] << 4);
}

uint8_t do_xor(uint8_t a, uint8_t b, uint8_t lookup_8x4[128])
{
	uint8_t ab = (a << 4) | b;
	uint8_t result;
	if (ab < 128) {
		result = lookup_8x4[ab] & 0x0F;
	}
	else
		result = (lookup_8x4[ab % 128] >> 4) & 0x0F;
	/* printf("do_xor: a=%x, b=%x, result=%x, a^b=%x\n", a, b, result, a^b); */
	return result;
}
uint8_t do_xor_hi_lo(uint8_t a, uint8_t b, uint8_t dual_lookup_8x4[][128])
{
	static const int hibits_lookup = 0;
	static const int lowbits_lookup = 1;
	uint8_t result = 0;
	int resultHi = do_xor(a>>4, b>>4, dual_lookup_8x4[hibits_lookup]);
	int resultLo = do_xor(a&0x0F, b&0x0F, dual_lookup_8x4[lowbits_lookup]);

	result = (resultHi << 4) | resultLo;
	/* printf("XORing HiLo: a=%x, b=%x, result=%x, a^b=%x\n", a, b, result, a^b); */
	return result;
}
void do_typeII(_4bit_strip32_t strips, uint8_t state[4][4],
		typeII_round_t round_typeIIs)
{
	int row, col, k;
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			/*  slice values must be re-arranged so that a simple xor can be done */
			/*  in the next step by means of lookup */
			for (k = 0; k < 4; ++k) {
				strips[k][col][row]
						= round_typeIIs[row][col][state[row][col]][k];
			}
		}
	}
}
void do_typeI(_4bit_strip128_t strips, uint8_t state[4][4], typeI_t typeIs)
{
	int row, col, k;
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			int index = row*4 + col;
			/* memcpy(strips[row][col], typeIs[row][col][state[row][col]], 16); */
			/* printf("strips[%d][%d] ", row, col); dump_hex("", strips[row][col], 16); */
			/*  slice values must be re-arranged so that a simple xor can be done */
			/*  in the next step by means of lookup */
			for (k = 0; k < 16; ++k) {
				strips[k/4][k%4][index]
						= typeIs[row][col][state[row][col]][k];
			}
		}
	}
}
void do_typeIV32(uint8_t state[4][4], _4bit_strip32_t strips,
		typeIV32_round_t round_typeIV32)
{
	int row, col, k;
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			k = 3;
			state[row][col] = strips[row][col][k];
			for (k = 2; k != -1; --k) {
				state[row][col] = do_xor_hi_lo(strips[row][col][k],
						state[row][col], round_typeIV32[row][col][k]);
			}
		}
	}
}
void do_typeIV128(uint8_t state[4][4], _4bit_strip128_t strips,
		typeIV128_t typeIV_Is)
{
	int row, col, k;
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			k = 15;
			state[row][col] = strips[row][col][k];
			for (k = 14; k != -1; --k) {
				state[row][col] = do_xor_hi_lo(strips[row][col][k],
						state[row][col], typeIV_Is[row][col][k]);
			}
		}
	}
}
void do_typeIV_II(uint8_t state[4][4], _4bit_strip32_t strips,
		typeIV_II_round_t round_typeIV_IIs)
{
	do_typeIV32(state, strips, round_typeIV_IIs);
}
void do_typeIV_III(uint8_t state[4][4], _4bit_strip32_t strips,
		typeIV_III_round_t round_typeIV_IIIs)
{
	do_typeIV32(state, strips, round_typeIV_IIIs);
}
void do_typeIII(_4bit_strip32_t strips, uint8_t state[4][4],
		typeIII_round_t round_typeIIIs)
{
	do_typeII(strips, state, round_typeIIIs);
}
void do_typeIA(_4bit_strip128_t strips, uint8_t state[4][4],
		typeIA_t typeIAs)
{
	do_typeI(strips, state, typeIAs);
}
void do_typeIB(_4bit_strip128_t strips, uint8_t state[4][4],
		typeIB_t typeIBs)
{
	do_typeI(strips, state, typeIBs);
}
void do_typeIV_IA(uint8_t state[4][4], _4bit_strip128_t strips,
		typeIV_IA_t typeIV_IAs)
{
	do_typeIV128(state, strips, typeIV_IAs);
}
void do_typeIV_IB(uint8_t state[4][4], _4bit_strip128_t strips,
		typeIV_IB_t typeIV_IBs)
{
	do_typeIV128(state, strips, typeIV_IBs);
}
void do_input(uint8_t state[4][4], const uint8_t input[KEY_SIZE],
		gf2matrix *linear_encoding, sboxes_8bit_t input_sboxes)
{
	int row, col, k;
	uint8_t v[16], w[16];
	/* reorder in a columnar manner */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			v[row*4 + col] = input[row + 4*col];
		}
	}
	mul_array_by_matrix_128x128(w, linear_encoding, v);
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			state[row][col] = sub_bytes_hi_lo(w[row*4 + col],
					input_sboxes[row][col]);
		}
	}
}
void do_output(uint8_t output[KEY_SIZE], uint8_t state[4][4],
		gf2matrix *linear_decoding, sboxes_8bit_t output_sboxes)
{
	int row, col;
	uint8_t v[16];
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			v[row*4 + col] = sub_bytes_hi_lo(state[row][col],
					output_sboxes[row][col]);
		}
	}
	mul_array_by_matrix_128x128(v, linear_decoding, v);
	/* reorder in a columnar manner */
	for (row = 0; row < 4; ++row) {
		for (col = 0; col < 4; ++col) {
			output[row + 4*col] = v[row* 4 + col];
		}
	}
}
