#ifndef INVERFC_H
#define INVERFC_H

#include "expr.h"

Expr* builtin_inverfc(Expr* res);   /* InverseErfc[s] */
void  inverfc_init(void);

#ifdef USE_MPFR
#include <mpfr.h>
/* erfcinv(s) for 0 < s < 2 into `out` (init2'd at `target` bits); robust for s
 * next to 0 or 2 and below the double range. Shared with InverseErf. */
void  inverfc_mpfr(mpfr_t out, const mpfr_t s, mpfr_prec_t target);
#endif

#endif /* INVERFC_H */
