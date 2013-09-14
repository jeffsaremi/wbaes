#include <stdio.h>
#include "aes.h"
#include "util.h"


uint8_t SBox[256] =
{
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5, 0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0, 0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC, 0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A, 0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0, 0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B, 0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85, 0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5, 0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17, 0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88, 0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C, 0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9, 0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6, 0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E, 0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94, 0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68, 0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};
uint8_t ISBox[256] =
{
    0x52, 0x09, 0x6A, 0xD5, 0x30, 0x36, 0xA5, 0x38, 0xBF, 0x40, 0xA3, 0x9E, 0x81, 0xF3, 0xD7, 0xFB,
    0x7C, 0xE3, 0x39, 0x82, 0x9B, 0x2F, 0xFF, 0x87, 0x34, 0x8E, 0x43, 0x44, 0xC4, 0xDE, 0xE9, 0xCB,
    0x54, 0x7B, 0x94, 0x32, 0xA6, 0xC2, 0x23, 0x3D, 0xEE, 0x4C, 0x95, 0x0B, 0x42, 0xFA, 0xC3, 0x4E,
    0x08, 0x2E, 0xA1, 0x66, 0x28, 0xD9, 0x24, 0xB2, 0x76, 0x5B, 0xA2, 0x49, 0x6D, 0x8B, 0xD1, 0x25,
    0x72, 0xF8, 0xF6, 0x64, 0x86, 0x68, 0x98, 0x16, 0xD4, 0xA4, 0x5C, 0xCC, 0x5D, 0x65, 0xB6, 0x92,
    0x6C, 0x70, 0x48, 0x50, 0xFD, 0xED, 0xB9, 0xDA, 0x5E, 0x15, 0x46, 0x57, 0xA7, 0x8D, 0x9D, 0x84,
    0x90, 0xD8, 0xAB, 0x00, 0x8C, 0xBC, 0xD3, 0x0A, 0xF7, 0xE4, 0x58, 0x05, 0xB8, 0xB3, 0x45, 0x06,
    0xD0, 0x2C, 0x1E, 0x8F, 0xCA, 0x3F, 0x0F, 0x02, 0xC1, 0xAF, 0xBD, 0x03, 0x01, 0x13, 0x8A, 0x6B,
    0x3A, 0x91, 0x11, 0x41, 0x4F, 0x67, 0xDC, 0xEA, 0x97, 0xF2, 0xCF, 0xCE, 0xF0, 0xB4, 0xE6, 0x73,
    0x96, 0xAC, 0x74, 0x22, 0xE7, 0xAD, 0x35, 0x85, 0xE2, 0xF9, 0x37, 0xE8, 0x1C, 0x75, 0xDF, 0x6E,
    0x47, 0xF1, 0x1A, 0x71, 0x1D, 0x29, 0xC5, 0x89, 0x6F, 0xB7, 0x62, 0x0E, 0xAA, 0x18, 0xBE, 0x1B,
    0xFC, 0x56, 0x3E, 0x4B, 0xC6, 0xD2, 0x79, 0x20, 0x9A, 0xDB, 0xC0, 0xFE, 0x78, 0xCD, 0x5A, 0xF4,
    0x1F, 0xDD, 0xA8, 0x33, 0x88, 0x07, 0xC7, 0x31, 0xB1, 0x12, 0x10, 0x59, 0x27, 0x80, 0xEC, 0x5F,
    0x60, 0x51, 0x7F, 0xA9, 0x19, 0xB5, 0x4A, 0x0D, 0x2D, 0xE5, 0x7A, 0x9F, 0x93, 0xC9, 0x9C, 0xEF,
    0xA0, 0xE0, 0x3B, 0x4D, 0xAE, 0x2A, 0xF5, 0xB0, 0xC8, 0xEB, 0xBB, 0x3C, 0x83, 0x53, 0x99, 0x61,
    0x17, 0x2B, 0x04, 0x7E, 0xBA, 0x77, 0xD6, 0x26, 0xE1, 0x69, 0x14, 0x63, 0x55, 0x21, 0x0C, 0x7D
};
static uint32_t Rcon[11] = {0x00000000,
	0x00000001, 0x00000002, 0x00000004, 0x00000008,
	0x00000010, 0x00000020, 0x00000040, 0x00000080,
	0x0000001B, 0x00000036
};

static uint8_t ModMul2[256];
static uint8_t ModMul3[256];
static uint8_t ModMul9[256];
static uint8_t ModMulB[256];
static uint8_t ModMulD[256];
static uint8_t ModMulE[256];
static int g_mod_mul = 0;
static int g_inv_mod_mul = 0;

static void gen_mod_mul()
{
	if(!g_mod_mul)
	{
		int i;
	    for(i = 0; i <= 256; ++i)
	    {
	        ModMul3[i] = (ModMul2[i] = xtime((uint8_t)i)) ^ (uint8_t)i;
	    }
	    /* dump_hex2("ModMul2: ", ModMul2, 256); */
	    /* dump_hex2("ModMul3: ", ModMul3, 256); */
	    g_mod_mul = 1;
	}
}
static void gen_inv_mod_mul()
{
	if(!g_inv_mod_mul)
	{
		int i;
	    for(i = 0; i <= 256; ++i)
	    {
	        uint8_t i1 = (uint8_t)i;
	        uint8_t i2 = xtime(i1);
	        uint8_t i4 = xtime(i2);
	        uint8_t i8 = xtime(i4);
	        ModMul9[i] = i8 ^ (uint8_t)i;
	        ModMulB[i] = i8 ^ i2 ^ i1;
	        ModMulD[i] = i8 ^ i4 ^ i1;
	        ModMulE[i] = i8 ^ i4 ^ i2;
	    }
	/*     dump_hex2("ModMul9: ", ModMul9, 256); */
	/*     dump_hex2("ModMulB: ", ModMulB, 256); */
	/*     dump_hex2("ModMulD: ", ModMulD, 256); */
	/*     dump_hex2("ModMulE: ", ModMulE, 256); */
	    g_inv_mod_mul = 1;
	}
}
void copy_word(uint8_t dest[4], uint8_t src[4])
{
    dest[0] = src[0];
    dest[1] = src[1];
    dest[2] = src[2];
    dest[3] = src[3];
}
void xor_words(uint8_t w1[4], uint8_t w2[4])
{
    w1[0] = w1[0] ^ w2[0];
    w1[1] = w1[1] ^ w2[1];
    w1[2] = w1[2] ^ w2[2];
    w1[3] = w1[3] ^ w2[3];
}
void add_round_key(uint8_t state[4][4], uint32_t key_schedule[4*(NR+1)],
		int round)
{
	int col;
    for(col = 0; col < 4; ++col)
    {
        uint8_t *key = (uint8_t*)&key_schedule[(round*4) + col];
        state[0][col] = state[0][col] ^ key[0];
        state[1][col] = state[1][col] ^ key[1];
        state[2][col] = state[2][col] ^ key[2];
        state[3][col] = state[3][col] ^ key[3];
    }
}
void mix_columns(uint8_t state[4][4])
{
	int col;
	if(!g_mod_mul)
		gen_mod_mul();
    for(col = 0; col < 4; ++col)
    {
        uint8_t temp[4];
        temp[0] = ModMul2[state[0][col]] ^ ModMul3[state[1][col]] ^ state[2][col] ^ state[3][col];
        temp[1] = state[0][col] ^ ModMul2[state[1][col]] ^ ModMul3[state[2][col]] ^ state[3][col];
        temp[2] = state[0][col] ^ state[1][col] ^ ModMul2[state[2][col]] ^ ModMul3[state[3][col]];
        temp[3] = ModMul3[state[0][col]] ^ state[1][col] ^ state[2][col] ^ ModMul2[state[3][col]];
        state[0][col] = temp[0];
        state[1][col] = temp[1];
        state[2][col] = temp[2];
        state[3][col] = temp[3];
    }
}
void mix_columns_inv(uint8_t state[4][4])
{
	int col;
	if(!g_inv_mod_mul)
		gen_inv_mod_mul();
	for(col = 0; col < 4; ++col)
    {
        uint8_t temp[4];
        temp[0] = ModMulE[state[0][col]] ^ ModMulB[state[1][col]] ^ ModMulD[state[2][col]] ^ ModMul9[state[3][col]];
        temp[1] = ModMul9[state[0][col]] ^ ModMulE[state[1][col]] ^ ModMulB[state[2][col]] ^ ModMulD[state[3][col]];
        temp[2] = ModMulD[state[0][col]] ^ ModMul9[state[1][col]] ^ ModMulE[state[2][col]] ^ ModMulB[state[3][col]];
        temp[3] = ModMulB[state[0][col]] ^ ModMulD[state[1][col]] ^ ModMul9[state[2][col]] ^ ModMulE[state[3][col]];
        state[0][col] = temp[0];
        state[1][col] = temp[1];
        state[2][col] = temp[2];
        state[3][col] = temp[3];
    }
}
void mix_columns_inv2(uint8_t key[4*4])
{
	int i;
	if(!g_inv_mod_mul)
		gen_inv_mod_mul();
    for(i = 0; i < 4*4; i+=4)
    {
        uint8_t temp[4];
        temp[0] = ModMulE[key[i]] ^ ModMulB[key[1+i]] ^ ModMulD[key[2+i]] ^ ModMul9[key[3+i]];
        temp[1] = ModMul9[key[i]] ^ ModMulE[key[1+i]] ^ ModMulB[key[2+i]] ^ ModMulD[key[3+i]];
        temp[2] = ModMulD[key[i]] ^ ModMul9[key[1+i]] ^ ModMulE[key[2+i]] ^ ModMulB[key[3+i]];
        temp[3] = ModMulB[key[i]] ^ ModMulD[key[1+i]] ^ ModMul9[key[2+i]] ^ ModMulE[key[3+i]];
        key[i] = temp[0];
        key[1+i] = temp[1];
        key[2+i] = temp[2];
        key[3+i] = temp[3];
    }
}
void rotate_word(uint8_t w[4], int i)
{
    uint8_t temp;
    switch(i)
    {
    case 1:
        {
            temp = w[0];
            w[0] = w[1];
            w[1] = w[2];
            w[2] = w[3];
            w[3] = temp;
            break;
        }
    case 2:
        {
            uint8_t temp2;
            temp = w[0]; temp2 = w[1];
            w[0] = w[2];
            w[1] = w[3];
            w[2] = temp;
            w[3] = temp2;
            break;
        }
    case 3:
        {
            temp = w[3];
            w[3] = w[2];
            w[2] = w[1];
            w[1] = w[0];
            w[0] = temp;
            break;
        }
    default:
        break;
    }
}
void shift_rows(uint8_t state[4][4])
{
    rotate_word((uint8_t*)&state[1][0], 1);
    rotate_word((uint8_t*)&state[2][0], 2);
    rotate_word((uint8_t*)&state[3][0], 3);
}
void inv_shift_rows(uint8_t state[4][4])
{
    rotate_word((uint8_t*)&state[1][0], 3);
    rotate_word((uint8_t*)&state[2][0], 2);
    rotate_word((uint8_t*)&state[3][0], 1);
}
void sub_word(uint8_t w[4], const uint8_t sbox[256])
{
	int i;
	for(i =0; i < 4; ++i)
		w[i] = sbox[(w[i])];
}
void sub_bytes(uint8_t state[4][4], const uint8_t sbox[256])
{
	int i;
    for(i = 0; i < 4; ++i)
        sub_word((uint8_t*)&state[i][0], sbox);
}
void inv_sub_bytes(uint8_t state[4][4], const uint8_t isbox[256])
{
	int i;
    for(i = 0; i < 4; ++i)
        sub_word((uint8_t*)&state[i][0], isbox);
}
void expand_key(uint8_t key[KEY_SIZE], const uint8_t sbox[256],
		uint32_t key_schedule[4*(NR+1)], int nk)
{
    int i = 0;
    while(i < nk)
    {
        copy_word((uint8_t*)&key_schedule[i], &key[4*i]);
        ++i;
    }
    i = nk;
    while(i < 4 * (NR+1))
    {
        uint8_t temp[4];
        copy_word(temp, (uint8_t*)&key_schedule[i-1]);
        /* printf("\n%02d:\t", i); */
        /* dump_word(temp, ""); */
        if(i % nk == 0)
        {
            rotate_word(temp, 1);
            /* dump_word(temp, ""); */
            sub_word(temp, sbox);
            /* dump_word(temp, ""); */
            /* dump_word((uint8_t*)&Rcon[i/nk], ""); */
            xor_words(temp, (uint8_t*)&Rcon[i/nk]);
            /* dump_word(temp, ""); */
        }
        else if((nk > 6) && (i % nk == 4))
        {
            /* printf("          "); */
            /* printf("          "); */
            /* printf("          "); */
            sub_word(temp, sbox);
            /* dump_word(temp, ""); */
        }
        else
        {
            /* printf("          "); */
            /* printf("          "); */
            /* printf("          "); */
            /* printf("          "); */
        }
        /* dump_word((uint8_t*)&key_schedule[i-nk], ""); */
        copy_word((uint8_t*)&key_schedule[i], (uint8_t*)&key_schedule[i-nk]);
        xor_words((uint8_t*)&key_schedule[i], temp);
        /* dump_word((uint8_t*)&key_schedule[i], ""); */
        ++i;
    }
    /* printf("\n"); */
}
/* void mix_expanded_key(uint32_t expanded_key[4*(NR+1)], */
/* 		uint32_t mixed_key_schedule[4*(NR+1)]) */
/* { */
/* 	int i, round; */
/*     for(i = 0; i < ((NR+1)*4); ++i) */
/*     { */
/*         mixed_key_schedule[i] = expanded_key[i]; */
/*     } */
/*     for(round = 1; round < NR; ++round) */
/*     { */
/*         mix_columns_inv2((uint8_t*)&mixed_key_schedule[round*4]); */
/*     } */
/* } */
void mix_expanded_key(uint32_t expanded_key[4*(NR+1)])
{
	int round;
    for(round = 1; round < NR; ++round)
    {
        mix_columns_inv2((uint8_t*)&expanded_key[round*4]);
    }
}
uint8_t xtime(uint8_t i)
{
    return (i << 1) ^ (( i & 0x80 ) ? 0x1B : 0x00 );
}
void input2state(uint8_t state[4][4], const uint8_t input[16])
{
	int i, j;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			state[i][j] = input[i + 4*j];
		}
	}
}
void state2output(uint8_t output[16], uint8_t state[4][4])
{
	int i, j;
	for (i = 0; i < 4; ++i)
	{
		for (j = 0; j < 4; ++j)
		{
			output[i + 4*j] = state[i][j];
		}
	}
}

void cipher(uint8_t in[4*4], uint8_t out[4*4], const uint8_t sbox[256],
		uint32_t key_schedule[4*(NR+1)])
{
	int round;
    uint8_t state[4][4];
    printf("%s\n", __FUNCTION__);


    dump_hex("Input: ", in, 4*4);
    input2state(state, in);
    dump_state("initial state: ", state);

    add_round_key(state, key_schedule, 0);

    for(round = 1; round < NR; ++round)
    {
/*         printf("round[%d]: \n", round); */
/*         dump_state("\t\tstart ", state); */
        sub_bytes(state, sbox);
/*         dump_state("\t\ts_box ", state); */
        shift_rows(state);
/*         dump_state("\t\ts_row ", state); */
        mix_columns(state);
/*         dump_state("\t\tm_col ", state); */
/*         dump_hex("\t\tk_sch ",(uint8_t*)&key_schedule[(round*4)], 4*4); */
        add_round_key(state, key_schedule, round);
    }
/*     printf("round[10]: \n"); */
/*     dump_state("\t\tstart ", state); */
    sub_bytes(state, sbox);
/*     dump_state("\t\ts_box ", state); */
    shift_rows(state);
/*     dump_state("\t\ts_row ", state); */
/*     dump_hex("\t\tk_sch ",(uint8_t*)&key_schedule[(10*4)], 4*4); */
    add_round_key(state, key_schedule, NR);

    dump_state("state after last round: ", state);

    state2output(out, state);
/*     dump_hex("Output: ", out, 4*4); */
}

void decipher(uint8_t in[4*4], uint8_t out[4*4], const uint8_t isbox[256],
		uint32_t key_schedule[4*(NR+1)])
{
	int round;
    uint8_t state[4][4];
    printf("%s\n", __FUNCTION__);

    /* dump_hex("Input: ", in, 4*4); */
    input2state(state, in);

    printf("round[10]: \n");
    dump_state("\t\tistart ", state);
/*     dump_hex("\t\tik_sch ",(uint8_t*)&key_schedule[10*4], 4*4); */
    add_round_key(state, key_schedule, 10);

    for(round = NR - 1; round >= 1; --round)
    {
        printf("round[%d]: \n", round);
        dump_state("\t\tistart ", state);
        inv_shift_rows(state);
/*         dump_state("\t\tis_row ", state); */
        inv_sub_bytes(state, isbox);
/*         dump_state("\t\tis_box ", state); */
/*         dump_hex("\t\tik_sch ",(uint8_t*)&key_schedule[(round*4)], 4*4); */
        add_round_key(state, key_schedule, round);
/*         dump_state("\t\tik_add ", state); */
        mix_columns_inv(state);
        /* dump_state("\t\tim_col ", state); */
    }
    printf("round[0]: \n");
    dump_state("\t\tistart ", state);
    inv_shift_rows(state);
/*     dump_state("\t\tis_row ", state); */
    inv_sub_bytes(state, isbox);
/*     dump_state("\t\tis_box ", state); */
/*     dump_hex("\t\tik_sch ",(uint8_t*)&key_schedule[0], 4*4); */
    add_round_key(state, key_schedule, 0);

    dump_state("state after last round: ", state);

    state2output(out, state);
/*     dump_hex("Output: ", out, 4*4); */
}

void eqv_decipher(uint8_t in[4*4], uint8_t out[4*4], const uint8_t isbox[256],
		uint32_t key_schedule[4*(NR+1)])
{
	int round;
    uint8_t state[4][4];
    printf("%s\n", __FUNCTION__);


    /* dump_hex("Input: ", in, 4*4); */
    input2state(state, in);

    printf("round[10]: \n");
    dump_state("\t\teistart ", state);
/*     dump_hex("\t\teik_sch ",(uint8_t*)&key_schedule[10*4], 4*4); */
    add_round_key(state, key_schedule, 10);

    for(round = NR - 1; round >= 1; --round)
    {
        printf("round[%d]: \n", round);
        dump_state("\t\teistart ", state);
        inv_sub_bytes(state, isbox);
/*         dump_state("\t\teis_box ", state); */
        inv_shift_rows(state);
/*         dump_state("\t\teis_row ", state); */
        mix_columns_inv(state);
/*         dump_state("\t\teim_col ", state); */
/*         dump_hex("\t\teik_sch ",(uint8_t*)&key_schedule[(round*4)], 4*4); */
        add_round_key(state, key_schedule, round);
    }
    printf("round[0]: \n");
    dump_state("\t\teistart ", state);
    inv_sub_bytes(state, isbox);
/*     dump_state("\t\teis_box ", state); */
    inv_shift_rows(state);
/*     dump_state("\t\teis_row ", state); */
/*     dump_hex("\t\tik_sch ",(uint8_t*)&key_schedule[0], 4*4); */
    add_round_key(state, key_schedule, 0);

    dump_state("state after last round: ", state);

    state2output(out, state);

/*     dump_hex("Output: ", out, 4*4); */
}
