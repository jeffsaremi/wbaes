#ifndef _WB_AES_GEN_H
#define _WB_AES_GEN_H
/**
 * \file wb_aes_gen.h
 * \brief WhiteBox AES key obfuscation functions
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */

#include <stdint.h>
#include "gf2matrix.h"
#include "wb_aes.h"

/** TBox Lookup table */
typedef uint8_t (tbox_t)[NR][4][4][256];
/** 32 4x4sboxes used in 4bit encoding/decoding */
typedef uint8_t (sboxes_128bit_t)[4][4][32][16];
/** 8 4x4sboxes used in 4bit encoding/decoding */
typedef uint8_t (sboxes_32bit_t)[4][4][8][16];
/** Mixing Bijections used along with TBoxes (in Type II) */
typedef gf2matrix *tbox_mixing_bijections_t[NR][4][4];

/**
 * \brief
 * substitute hi and low 4-bits in x using the given SBOXes
 * for a given 4-byte array
 */
void sub_bytes_hi_lo32(uint8_t v[4], uint8_t _sbox[][16]);
/**
 * \brief
 * substitute hi and low 4-bits in x using the given SBOXes
 * for a given 16-byte array
 */
void sub_bytes_hi_lo128(uint8_t v[16], uint8_t _sbox[][16]);
/**
 * \brief
 * TBoxes are used in Type II construction
 * _tbox is a combination of original AES _sbox with the expanded round key applied
 * A mixing bijection is applied (also see Type III for the inverse of this m.b.)
 * Note:
 * Shift Rows step is moved up to before SubBytes (which is OK) and to prior to
 * Add_roundKey step. This means the key that we select should be the key prior to
 * the shift. In other terms we need to apply an inverse shift to the state byte
 * index to find the index for the corresponding key byte
 */
int make_tbox(tbox_t tbox,
		uint8_t sbox[256],
		uint32_t expanded_key[(NR+1)*4],
		tbox_mixing_bijections_t tbox_mixing_bijection);
/**
 * \brief
 * Same as make_tbox but for decryption
 */
int make_inv_tbox(tbox_t tbox,
		uint8_t sbox[256],
		uint32_t expanded_key[(NR+1)*4],
		tbox_mixing_bijections_t tbox_mixing_bijection);
/**
 * \brief
 * Makes Type IA lookup tables corresponding to each byte in the state
 * \param typeIA	The Type IA lookup table
 * \param first_inv_tbox_mixing_bijection The inverse of the first TBox Mixing Bijection
 * which was used in Type II in the first round
 * \param initial_decoding The inverse of a corresponding 128x128 diffusion matrix.
 * The encoding side will be embedded in the client.
 * \param decoding_sbox Inverses of non-linear sboxes at the input to this Type.
 * The encoding side is embedded in the caller.
 * \param encoding_sbox Non-linear sboxes used to encode the output from this round.
 * The corresponding inverses will be used at the input to the next step (TypeIV)
 */
int make_typeIA(typeIA_t typeIA,
		gf2matrix *first_inv_tbox_mixing_bijection[4][4],
		gf2matrix *initial_decoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox);
/**
 * \brief
 * Makes Type IB lookup tables corresponding to each byte in the state
 * \param typeIB	Type IB lookup table
 * \param last_round_tbox The TBoxes corresponding to the last two rounds (9+10)
 * which was used in Type II in the first round
 * \param final_encoding A 128x128 diffusion matrix.
 * The deccoding side will be embedded in the client.
 * \param decoding_sbox Inverses of non-linear sboxes at the input to this Type.
 * The encoding side is embedded in the caller.
 * \param encoding_sbox Non-linear sboxes used to encode the output from this round.
 * The corresponding inverses will be used at the input to the next step (TypeIV)
 */
int make_typeIB(typeIB_t typeIB,
		uint8_t last_round_tbox[4][4][256],
		gf2matrix *final_encoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox);
/**
 * \brief
 * Same as make_typeIB but for decryption
 */
int make_inv_typeIB(typeIB_t typeIB,
		uint8_t last_round_tbox[4][4][256],
		gf2matrix *final_encoding,
		sboxes_8bit_t decoding_sbox,
		sboxes_128bit_t encoding_sbox);
/**
 * \brief
 * Type II tables are a combination of TBoxes, original AES Mix _columns,
 * and a mixing bijection on top of them
 */
int make_typeII(typeII_t typeII,
		tbox_t tbox,
		gf2matrix *mix_columns_mixing_bijection,
		sboxes_8bit_t decoding_sbox[NR-1],
		sboxes_32bit_t encoding_sbox[NR-1]);
/**
 * \brief
 * Same as make_typeII but for decryption
 */
int make_inv_typeII(typeII_t typeII,
		tbox_t tbox,
		gf2matrix *mix_columns_mixing_bijection,
		sboxes_8bit_t decoding_sbox[NR-1],
		sboxes_32bit_t encoding_sbox[NR-1]);
/**
 * \brief
 * Type III tables are to reverse the diffusion introduced in Type II
 * \param typeIII The lookup tables for this round (output)
 * \param inv_mix_columns_mixing_bijection	The inverse of the additional mixing bijetion used in Type II
 * \param inv_tbox_mixing_bijections	The inverses of the mising bijections used in previous Type II
 * \param decoding_sbox Inverses of non-linear sboxes at the input to this Type.
 * The encoding sides were used in the previous Type IV step.
 * \param encoding_sbox Non-linear sboxes used to encode the output from this round.
 * The correspnding decoding sboxes will be used in the next Type IV step.
 */
int make_typeIII(typeIII_t typeIII,
		gf2matrix *inv_mix_columns_mixing_bijection,
		tbox_mixing_bijections_t inv_tbox_mixing_bijections,
		sboxes_8bit_t decoding_sbox[NR-1],
		sboxes_32bit_t encoding_sbox[NR-1]);
/**
 * \brief
 * A (2^8)x4 bit lookup or a 256x4bits is coded like a 128 byte array.
 * The low 4 bits in each array item represents the result for the first 128.
 * The hi 4 bits are for the next 128 (128-255).
 * To lookup the Xor of "a" (4bits) and "b" (4bits) the following must be done:
 * 1. create an index by combining a and b: (a<<4)|b. This is an 8bit index.
 * 2. If the index is < 128 then the result of the lookup would be the
 * Table[index] & 0x0F (the lo 4 bits)
 * 3. If the index is >= 128 then the result of the lookup would be the hi 4bits:
 * (Table[index%128] >> 4) & 0x0F
 */
int make_typeIV128(typeIV128_t typeIV128,
		sboxes_128bit_t decoding_sbox,
		sboxes_8bit_t encoding_sbox,
		sboxes_8bit_t encoding_sbox_inv);
/**
 * \brief
 * \see make_typeIV128
 */
int make_typeIV_IA(typeIV128_t typeIV_IA,
		sboxes_128bit_t decoding_sbox,
		sboxes_8bit_t encoding_sbox,
		sboxes_8bit_t encoding_sbox_inv);
/**
 * \brief
 * \see make_typeIV128
 */
int make_typeIV_IB(typeIV128_t typeIV_IB,
		sboxes_128bit_t decoding_sbox,
		sboxes_8bit_t encoding_sbox,
		sboxes_8bit_t encoding_sbox_inv);
/**
 * \brief
 * make_typeIV32 prepares lookup tables for TypeIV transformations
 * after TypeII and after TypeIII.
 * A (2^8)x4 bit lookup or a 256bytesx4bits is coded like a 128 byte array.
 * The low 4 bits in each array item represents the result for the first 128.
 * The hi 4 bits are for the next 128 (128-255).
 * To lookup the Xor of a (4bits) and b (4bits) the following must be done:
 * 1. create an index by combining a and b: (a<<4)|b. This is an 8bit index.
 * 2. If the index is < 128 then the result of the lookup would be the
 * Table[index] & 0x0F (the lo 4 bits)
 * 3. If the index is >= 128 then the result of the lookup would be the hi 4bits:
 * (Table[index%128] >> 4) & 0x0F
 */
int make_rounds_typeIV(typeIV32_t typeIV32,
		sboxes_32bit_t decoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox_inv[NR - 1]);
/**
 * \see make_rounds_typeIV
 */
int make_typeIV_II(typeIV32_t typeIV_II,
		sboxes_32bit_t decoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox_inv[NR - 1]);
/**
 * \see make_rounds_typeIV
 */
int make_typeIV_III(typeIV32_t typeIV_III,
		sboxes_32bit_t decoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox[NR - 1],
		sboxes_8bit_t encoding_sbox_inv[NR - 1]);
/**
 * \brief
 * Slice calculation using the original AES Mix-Columns matrix
 */
void calc_mixing_slices(uint8_t mc[4],
		uint8_t tbox_val,
		int slice_num);
/**
 * \brief
 * Inverse slice calculation using the original AES Inverse Mix-Columns matrix
 */
void calc_inv_mixing_slices(uint8_t mc[4],
		uint8_t tbox_val,
		int slice_num);

/*
 * ******** Helper Functions to create multi-dimensional arrays
 */
/**
 * \brief
 * populate the gives arrays with 8 4x4 Sboxes for the given number of rounds
 */
int make_rounds_sbox_pair_32(sboxes_32bit_t *encoding_sboxes,
		sboxes_32bit_t *decoding_sboxes, int rounds);
/**
 * \brief
 * populate the gives arrays with 8 4x4 Sboxes for the given number of rounds
 */
int make_rounds_sbox_pair_8(sboxes_8bit_t *encoding_sboxes,
		sboxes_8bit_t *decoding_sboxes, int rounds);
/**
 * \brief
 * populate the gives arrays with 32 4x4 Sboxes
 */
int make_sbox_pair_128(sboxes_128bit_t encoding_sboxes,
		sboxes_128bit_t decoding_sboxes);
/**
 * \brief
 * populate the gives arrays with 8 4x4 Sboxes
 */
int make_sbox_pair_32(sboxes_32bit_t encoding_sboxes,
		sboxes_32bit_t decoding_sboxes);
/**
 * \brief
 * populate the gives arrays with 2 4x4 Sboxes
 */
int make_sbox_pair_8(sboxes_8bit_t encoding_sboxes,
		sboxes_8bit_t decoding_sboxes);

/**
 * \brief
 * Utility to populate mixing bijections (and their inverses) used in TBox generation
 */
int make_tbox_mixing_bijections(tbox_mixing_bijections_t tbox_mbs,
		tbox_mixing_bijections_t tbox_mb_invs);
/**
 * \brief
 * Frees the matrices allocated in make_tbox_mixing_bijections()
 */
void free_tbox_mixing_bijections(tbox_mixing_bijections_t tbox_mbs);

#endif /*  _WB_AES_GEN_H */
