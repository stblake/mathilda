#ifndef PRODUCT_H
#define PRODUCT_H

#include "expr.h"

/*
 * Product[f, {i, imax}]                evaluate prod_{i=1}^{imax} f
 * Product[f, {i, imin, imax}]          prod_{i=imin}^{imax} f
 * Product[f, {i, imin, imax, di}]      step di
 * Product[f, {i, {i1, i2, ...}}]       successive values
 * Product[f, {i,...}, {j,...}, ...]    multiple (nested) products
 * Product[f, i]                        indefinite product (anti-quotient)
 *
 * The multiplicative analogue of Sum (src/sum/).  Implemented as a dispatcher
 * (this file) that handles iteration and a Method cascade over the
 * context-qualified sub-algorithms Product`Telescoping, Product`Rational,
 * Product`Geometric, Product`QProduct, ... -- mirroring how Sum dispatches to
 * Sum`Polynomial etc. and Integrate to Integrate`BronsteinRational.
 */
Expr* builtin_product(Expr* res);

void product_init(void);

#endif
