#ifndef M4RI_G2MATRIX_TYPEDEF_H_
#define M4RI_G2MATRIX_TYPEDEF_H_
/**
 * \file m4ri_gf2matrix_typedef.h
 * \brief GF(2) Matrix Typedefs.
 * This links GF2Matrix definition with that of M4RI's
 *
 * \author Jeff Saremi <jeffsaremi@yahoo.com>
 */
#define rotate_word m4ri_rotate_word
#include "m4ri/m4ri.h"
#undef rotate_word

/** make gf2matrix the same as mzd_t in m4ri */
typedef mzd_t gf2matrix;

#endif /* M4RI_G2MATRIX_TYPEDEF_H_ */
