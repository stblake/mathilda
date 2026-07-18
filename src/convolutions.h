#ifndef MATHILDA_CONVOLUTIONS_H
#define MATHILDA_CONVOLUTIONS_H

/* ---------------------------------------------------------------------------
 * ListConvolve / ListCorrelate — discrete convolution and correlation.
 *
 * ListConvolve[ker, list] forms Sum_r ker[r] list[s-r]; ListCorrelate forms
 * Sum_r ker[r] list[s+r]. Both support the full Wolfram-Language argument set:
 * cyclic overhangs {kL, kR} (or a single index k), constant / cyclic / empty
 * padding, generalized g / h in place of Times / Plus, and multidimensional
 * kernels and data. Data may be symbolic, exact, machine-precision, or
 * arbitrary-precision (MPFR).
 *
 * A fully general direct engine is the correctness backbone; for large numeric
 * inputs with the default Times/Plus a separable FFT fast path (FFTW for
 * machine precision, the MPFR FFT for arbitrary precision) is used instead, in
 * both 1-D and n-D. The two builtins share one engine, parameterised by mode.
 * -------------------------------------------------------------------------- */

#include "expr.h"

typedef enum { CONV_MODE_CONVOLVE, CONV_MODE_CORRELATE } ConvMode;

/* Shared engine entry point. `res` is the borrowed Head[args...] expression. */
Expr* conv_engine(Expr* res, ConvMode mode);

Expr* builtin_list_convolve(Expr* res);
void  convolutions_init(void);

#endif /* MATHILDA_CONVOLUTIONS_H */
