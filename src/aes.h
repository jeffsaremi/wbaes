#ifndef _AES_H_
#define _AES_H_
/**
 * \file aes.h
 * \brief Text-book Implementation of AES
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */

#include <stdint.h>

/** Number of rounds (set to 10 for AES128) */
#define NR 10
/** AES key size (set to 16 for AES128) */
#define KEY_SIZE 16

/** The original forward SBOX of AES */
extern uint8_t SBox[256];
/** The original inverse SBOX of AES */
extern uint8_t ISBox[256];
/** verbatim from AES standard */
void copy_word(uint8_t dest[4], uint8_t src[4]);
/** verbatim from AES standard */
void xor_words(uint8_t w1[4], uint8_t w2[4]);
/** verbatim from AES standard */
void add_round_key(uint8_t state[4][4], uint32_t key_schedule[4*(NR+1)],
		int round);
/** verbatim from AES standard */
void mix_columns(uint8_t state[4][4]);
/** verbatim from AES standard */
void mix_columns_inv(uint8_t state[4][4]);
/** verbatim from AES standard */
void mix_columns_inv2(uint8_t key[4*4]);
/** verbatim from AES standard */
void rotate_word(uint8_t w[4], int i);
/** verbatim from AES standard */
void shift_rows(uint8_t state[4][4]);
/** verbatim from AES standard */
void inv_shift_rows(uint8_t state[4][4]);
/** verbatim from AES standard */
void sub_word(uint8_t w[4], const uint8_t sbox[256]);
/** verbatim from AES standard */
void sub_bytes(uint8_t state[4][4], const uint8_t sbox[256]);
/** verbatim from AES standard */
void inv_sub_bytes(uint8_t state[4][4], const uint8_t isbox[256]);
/** verbatim from AES standard */
uint8_t xtime(uint8_t i);
/** Utility to copy 16 bytes of input to the state matrix */
void input2state(uint8_t state[4][4], const uint8_t input[16]);
/** Utility to copy state bytes to the 16 bytes of output */
void state2output(uint8_t output[16], uint8_t state[4][4]);
/** verbatim from AES standard */
void expand_key(uint8_t key[KEY_SIZE], const uint8_t sbox[256],
		uint32_t key_schedule[4*(NR+1)], int nk);
/** verbatim from AES standard */
void mix_expanded_key(uint32_t expanded_key[4*(NR+1)]);
/** cipher implementation as per AES standard */
void cipher(uint8_t in[4*4], uint8_t out[4*4], const uint8_t sbox[256],
		uint32_t key_schedule[4*(NR+1)]);
/** Inverse cipher implementation as per AES standard */
void decipher(uint8_t in[4*4], uint8_t out[4*4], const uint8_t isbox[256],
		uint32_t key_schedule[4*(NR+1)]);
/** Equivalent Inverse cipher implementation as per AES standard */
void eqv_decipher(uint8_t in[4*4], uint8_t out[4*4], const uint8_t isbox[256],
		uint32_t key_schedule[4*(NR+1)]);

#endif /* _AES_H_ */
