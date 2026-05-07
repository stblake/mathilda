#ifndef ARITHMETIC_H
#define ARITHMETIC_H

#include "expr.h"

// GCD helper
int64_t gcd(int64_t a, int64_t b);

// LCM helper
int64_t lcm(int64_t a, int64_t b);

// Create a simplified Rational or Integer
Expr* make_rational(int64_t n, int64_t d);

// Check if expression is a Rational[...]
bool is_rational(Expr* e, int64_t* n, int64_t* d);

// Check if expression is a Complex[re, im]
bool is_complex(Expr* e, Expr** re, Expr** im);

// Create a Complex expression
Expr* make_complex(Expr* re, Expr* im);

// Infinity / Indeterminate predicates
bool is_infinity_sym(Expr* e);
bool is_complex_infinity_sym(Expr* e);
bool is_indeterminate_sym(Expr* e);
// True if e == Times[c, Infinity] with c a negative numeric.
bool is_neg_infinity_form(Expr* e);
// Sign of a numeric expression (Integer, Real, BigInt, Rational[n,d]).
// Returns +1, -1, 0 for known signs; 0 for non-numeric or ambiguous.
int  expr_numeric_sign(Expr* e);

/* expr_is_superficially_negative:
 * True if `e` can be cheaply rewritten with its sign flipped. Handles:
 *   - negative Integer/Real/Rational[-n, d]
 *   - Complex[re, im] where (sign(re) < 0) or (re == 0 and sign(im) < 0)
 *   - Times[c, ...] whose leading factor c is itself superficially negative
 * Used by trig/hyperbolic sign-extraction to pull out the sign via
 * f[-x] = +/- f[x] for even/odd f. `e` is never modified. */
bool expr_is_superficially_negative(Expr* e);

/* Re-entrant mute for arithmetic warnings (Power::infy, Infinity::indet, ...).
 * Internal probes that evaluate candidate sub-expressions — e.g. Limit — bump
 * this counter while they poke at divergent forms, then decrement on exit.
 * When the counter is non-zero, arithmetic modules skip their stderr messages.
 * The messages are informational only; the return value (ComplexInfinity,
 * Indeterminate) is unchanged. */
extern int g_arith_warnings_muted;
static inline void arith_warnings_mute_push(void)   { g_arith_warnings_muted++; }
static inline void arith_warnings_mute_pop(void)    { if (g_arith_warnings_muted > 0) g_arith_warnings_muted--; }
static inline int  arith_warnings_muted(void)       { return g_arith_warnings_muted; }

Expr* builtin_divide(Expr* res);
Expr* builtin_subtract(Expr* res);
Expr* builtin_complex(Expr* res);
Expr* builtin_rational(Expr* res);
Expr* builtin_gcd(Expr* res);
Expr* builtin_add(Expr* res);
Expr* builtin_lcm(Expr* res);
Expr* builtin_powermod(Expr* res);
Expr* builtin_factorial(Expr* res);
Expr* builtin_factorial2(Expr* res);
Expr* builtin_binomial(Expr* res);

#endif
