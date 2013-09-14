#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include <stdint.h>
#include "aes.h"
#include "util.h"
#include "test.h"
#include "test_util.h"

extern uint8_t SBox[256];
extern uint8_t ISBox[256];

void test_cipher_decipher()
{
	uint32_t key_schedule[4*(NR+1)] = {0};
	uint8_t key[KEY_SIZE] = {
		0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6,
		0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
	uint8_t in[4*4] = {
		0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
		0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
	uint8_t out[4*4], decipher_out[4*4];

	expand_key(key, SBox, key_schedule, 4);
	cipher(in, out, SBox, key_schedule);
	decipher(out, decipher_out, ISBox, key_schedule);

	ASSERT(memcmp(in, decipher_out, 16) == 0);
}
void test_cipher_decipher_eq()
{
	uint32_t key_schedule[4*(NR+1)] = {0};
	uint8_t key[KEY_SIZE] = {
		0x2b, 0x7e, 0x15, 0x16, 0x28, 0xae, 0xd2, 0xa6,
		0xab, 0xf7, 0x15, 0x88, 0x09, 0xcf, 0x4f, 0x3c};
	uint8_t in[4*4] = {
		0x00, 0x11, 0x22, 0x33, 0x44, 0x55, 0x66, 0x77,
		0x88, 0x99, 0xaa, 0xbb, 0xcc, 0xdd, 0xee, 0xff};
	uint8_t out[4*4], eq_decipher_out[16];

	expand_key(key, SBox, key_schedule, 4);
	cipher(in, out, SBox, key_schedule);
	mix_expanded_key(key_schedule);
	eqv_decipher(out, eq_decipher_out, ISBox, key_schedule);

	ASSERT(memcmp(in, eq_decipher_out, 16) == 0);
}
int main()
{
	TEST(test_cipher_decipher);
	TEST(test_cipher_decipher_eq);
	return 0;
}
