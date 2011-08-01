#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "util.h"

int get_random(int min, int max)
{
    if(max < min)
    {
    	int t = max;
    	max = min;
    	min = t;
    }
    if(max == 0)
        return 0;

    return min + (rand() % ((max-min)+1));
}
void dump_4bit_strip32(const char *label, uint8_t _4bitStrip[4][4][4])
{
	int i, j, k;
    if(label)
        printf("%s", label);
    printf("\n");
    for(i = 0; i < 4; ++i)
    {
    	printf("row %d: ", i);
        for(j = 0; j < 4; ++j)
        {
        	printf(" ");
        	for(k = 0; k < 4; ++k)
        	{
        		printf("%02x", _4bitStrip[i][j][k]);
        	}
        }
    	printf("\n");
    }
    printf("\n");
}
void dump_4bit_strip128(const char *label, uint8_t _4bitStrip[4][4][16])
{
	int i, j, k;
    if(label)
        printf("%s", label);
    printf("\n");
    for(i = 0; i < 4; ++i)
    {
    	printf("row %d: ", i);
        for(j = 0; j < 4; ++j)
        {
        	printf(" ");
        	for(k = 0; k < 16; ++k)
        	{
        		printf("%02x", _4bitStrip[i][j][k]);
        	}
        }
    	printf("\n");
    }
    printf("\n");
}
void dump_hex(const char *label, const uint8_t *buf, int size)
{
	int i;
    if(label)
        printf("%s", label);
    for(i = 0; i < size; ++i)
        printf("%02x", buf[i]);
    printf("\n");
}
void dump_hex2(const char *label, const uint8_t *buf, int size)
{
	int i, j;
    int indent = label ? strlen(label) : 0;
    if(label)
        printf("%s", label);
    for(i = 0; i < size; ++i)
    {
        if( (i%16) == 0 )
        {
            printf("\n");
            for(j = 0; j < indent; ++j)
                printf(" ");
        }
        printf("%02x", buf[i]);
    }
    printf("\n");
}
void dump_state(const char *label, uint8_t state[4][4])
{
	int i, j;
    if(label)
        printf("%s", label);
    for(i = 0; i < 4; ++i)
    {
        for(j = 0; j < 4; ++j)
            printf("%02x", state[i][j]);
        printf(" ");
    }
    printf("\n");
}
void dump_word(const uint8_t w[4], const char *label)
{
	int i;
    if(label)
        printf("%s", label);
    for(i = 0; i < 4; ++i)
        printf("%02x", w[i]);
    printf("  ");
}
void dump_decimal(const char *label, const uint8_t *buf, int size)
{
	int i;
    if(label)
        printf("%s", label);
    for(i = 0; i < size; ++i)
        printf("%d ", buf[i]);
    printf("\n");
}
