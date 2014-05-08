#include <stdint.h>

/**
 * Helper macros to compute cost of our computations.
 *
 * Usage:
 *   COST_MEASURE
 *     If defined, counters are activated.
 *     Otherwise, no code is added to source.
 *     NOTE: this must be defined before the header is included
 *
 *   COST_VARIABLES_HERE
 *     Declare this once for your project, e.g. in the main file.
 *
 *   COST_INC_yyy(x)
 *     Increases counter 'yyy' by value 'x'.
 *
 *   COST_yyy
 *     Value of counter 'yyy'.
 *
 *   Valid counters are: ADD, MUL, DIV, POW, ABS.
 */

#define COST_T   uint64_t
#define PRI_COST PRId64

#ifndef COST_MEASURE

#define COST_VARIABLES_HERE
#define COST_MODEL_RESET

#define COST_ADD ((const COST_T) 0)
#define COST_MUL ((const COST_T) 0)
#define COST_DIV ((const COST_T) 0)
#define COST_POW ((const COST_T) 0)
#define COST_ABS ((const COST_T) 0)

#define COST_INC_ADD(x)
#define COST_INC_MUL(x)
#define COST_INC_DIV(x)
#define COST_INC_POW(x)
#define COST_INC_ABS(x)

#else

extern COST_T __cost_add;
extern COST_T __cost_mul;
extern COST_T __cost_div;
extern COST_T __cost_pow;
extern COST_T __cost_abs;

#define COST_VARIABLES_HERE  \
    COST_T __cost_add = 0;   \
    COST_T __cost_mul = 0;   \
    COST_T __cost_div = 0;   \
    COST_T __cost_pow = 0;   \
    COST_T __cost_abs = 0;

#define COST_MODEL_RESET \
    __cost_add = 0;      \
    __cost_mul = 0;      \
    __cost_div = 0;      \
    __cost_pow = 0;      \
    __cost_abs = 0;

#define COST_ADD ((const COST_T) __cost_add)
#define COST_MUL ((const COST_T) __cost_mul)
#define COST_DIV ((const COST_T) __cost_div)
#define COST_POW ((const COST_T) __cost_pow)
#define COST_ABS ((const COST_T) __cost_abs)

#define COST_INC_ADD(x) {__cost_add += (x);}
#define COST_INC_MUL(x) {__cost_mul += (x);}
#define COST_INC_DIV(x) {__cost_div += (x);}
#define COST_INC_POW(x) {__cost_pow += (x);}
#define COST_INC_ABS(x) {__cost_abs += (x);}

#endif
