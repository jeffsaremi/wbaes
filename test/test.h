#ifndef TEST_H_
#define TEST_H_

#include <stdio.h>
#include <string.h>

#define ASSERT(x) if(!x) { \
	fprintf(stdout, "\n\n\t*** Failed: %s, %s, %s:%d\n\n", \
			#x, __FUNCTION__, __FILE__, __LINE__); \
	exit(1); }

#define ASSERT_EQUAL(a, b) if((a) != (b)) { \
	fprintf(stdout, "\n\n\tFailed: expected=%d(0x%x), got=%d(0x%x)" \
			" inside %s, %s:%d\n\n", \
			(a), (a), (b), (b), __FUNCTION__, __FILE__, __LINE__); \
	exit(1); }

#define TEST(name) { \
    printf("******************************************\n"); \
    printf("**** Test: " #name "\n");     \
    printf("****\n"); \
    name();}


#endif /* TEST_H_ */
