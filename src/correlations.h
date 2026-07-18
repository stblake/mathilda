#ifndef MATHILDA_CORRELATIONS_H
#define MATHILDA_CORRELATIONS_H

/* ---------------------------------------------------------------------------
 * ListCorrelate — discrete correlation. A thin front end over the shared
 * convolution/correlation engine in convolutions.c (see convolutions.h).
 * -------------------------------------------------------------------------- */

#include "expr.h"

Expr* builtin_list_correlate(Expr* res);
void  correlations_init(void);

#endif /* MATHILDA_CORRELATIONS_H */
