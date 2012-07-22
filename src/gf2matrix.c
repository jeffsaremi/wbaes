#include <stdio.h>
#include <stdlib.h>
#include "gf2matrix.h"
#include "util.h"

int is_invertible(const gf2matrix *m)
{
	int rc = 0;
	gf2matrix *inv = invert_matrix(NULL, m);
	if (inv) {
		free_matrix(inv);
		rc = 1;
	}

	return rc;
}
int comp_regions(const gf2matrix *a, int start_row_a, int start_col_a,
		const gf2matrix *b, int start_row_b, int start_col_b, int row_count,
		int col_count)
{
	int i, j;
	for (i = 0; i < row_count; ++i) {
		for (j = 0; j < col_count; ++j) {
			if (get_bit_at(a, start_row_a + i, start_col_a + j) != get_bit_at(
					b, start_row_b + i, start_col_b + j))
				return -1;
		}
	}
	return 0;
}

int copy_matrix(gf2matrix *dest, const gf2matrix *src)
{
	if (!dest || !src)
		return -1;
	if ((get_rows(src) != get_rows(dest)) || (get_cols(src) != get_rows(dest)))
		return -2;
	return copy_matrix_to_offset(dest, src, 0, 0);
}
