/*
 * wavelets.h -- HaarWavelet[], DaubechiesWavelet[n] and their scaling /
 * wavelet functions WaveletPhi, WaveletPsi, plus WaveletFilterCoefficients
 * (src/wavelets.c).  See docs/spec/builtins/special-functions.md.
 */
#ifndef MATHILDA_WAVELETS_H
#define MATHILDA_WAVELETS_H

#include "expr.h"

void  wavelets_init(void);
Expr* builtin_waveletphi(Expr* res);
Expr* builtin_waveletpsi(Expr* res);
Expr* builtin_waveletfiltercoefficients(Expr* res);
Expr* builtin_daubechieswavelet(Expr* res);

#endif /* MATHILDA_WAVELETS_H */
