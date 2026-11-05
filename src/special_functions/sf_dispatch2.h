#ifndef MATHILDA_SF_DISPATCH2_H
#define MATHILDA_SF_DISPATCH2_H

/* The numeric dispatch shared by the two-argument special functions f[nu, z]
 * whose machine kernel is a raw-double helper and whose arbitrary-precision /
 * complex path is an arb bridge (StruveL, AngerJ, WeberE, the Kelvin family).
 * It is StruveH's dispatch, factored out:
 *
 *   no inexact argument, or a non-numeric one   -> *handled = false
 *   target precision > 53 bits                  -> arb at that precision
 *   machine precision, both arguments real      -> machine(which, nu, z) -- the
 *                                                  SAME helper the NDArray kernel
 *                                                  and Compile[] call, so the
 *                                                  surfaces agree bit for bit
 *   machine kernel declined / complex argument  -> arb at 53 bits, demoted to
 *                                                  machine numbers
 *
 * Returns the value, or NULL (with *handled = true) to leave the call
 * unevaluated (arb declined: a pole, or FLINT absent). */

#include "expr.h"
#include <stdbool.h>

typedef bool  (*sf2_machine_fn)(int which, double nu, double z, double* out);
typedef Expr* (*sf2_arb_fn)(const Expr* nu, const Expr* z, int which);

Expr* sf2_numeric(const Expr* nu, const Expr* z, int which,
                  sf2_machine_fn machine, sf2_arb_fn arb, bool* handled);

/* True when e is an exact integer (Integer or bignum) -- the "exact integer
 * order" test of the closed-form reductions. */
bool sf2_exact_integer(const Expr* e);

/* True when e is the exact Integer 0. */
bool sf2_exact_zero(const Expr* e);

#endif /* MATHILDA_SF_DISPATCH2_H */
