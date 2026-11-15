/*
 * ff_fnx.h -- the field K = Q(v_1, ..., v_r)[t]/(mu(t)) as a FLINT generic
 * ring (gr_ctx), part of the ffield library (src/poly/ffield/, Phase 4.5 of
 * the Bronstein completion campaign: efficient ideal and field arithmetic in
 * C, reusable outside the integrators).
 *
 * mu is a monic irreducible polynomial over Q of degree n >= 1; n = 1 (mu = t)
 * is the plain rational function field Q(v_1..v_r).  Since Q is algebraically
 * closed in Q(v), mu stays irreducible over Q(v) and K is a field: K = Q(theta)(v)
 * with theta a root of mu.  This is Trager's fieldFn with an extension.
 *
 * An element is the dense coefficient vector (a_0, ..., a_(n-1)) of
 * a_0 + a_1 t + ... + a_(n-1) t^(n-1), each a_j an fmpz_mpoly_q in canonical
 * form (numerator and denominator coprime, denominator with positive leading
 * coefficient), so equality is coefficientwise and exact.
 *
 * Every gr_poly / gr_mat algorithm of FLINT then runs over K unchanged:
 * Hermite forms over K[x], gcds, determinants (src/poly/flint_polymat.c).
 *
 * Without FLINT the header declares nothing usable; callers are compiled under
 * USE_FLINT only.
 */
#ifndef FF_FNX_H
#define FF_FNX_H

#ifdef USE_FLINT
#include <flint/flint.h>
#include <flint/fmpq.h>
#include <flint/fmpq_poly.h>
#include <flint/fmpz_mpoly.h>
#include <flint/fmpz_mpoly_q.h>
#include <flint/gr.h>

typedef struct {
    fmpz_mpoly_ctx_t mctx;   /* Q[v_1..v_r] (r >= 1; a dummy variable when r = 0) */
    slong nvars;             /* r */
    slong n;                 /* deg mu */
    fmpq* mu;                /* mu_0 .. mu_n, mu_n = 1 */
} ff_fnx_data;

#define FF_FNX_DATA(ctx) ((ff_fnx_data*) GR_CTX_DATA_AS_PTR(ctx))

/* Initialise ctx as Q(v_1..v_nvars)[t]/(mu).  mu monic of degree >= 1 over Q
 * (irreducibility is the caller's responsibility; NULL means mu = t, n = 1).
 * Returns GR_SUCCESS, or GR_DOMAIN for a mu that is not monic. */
int ff_fnx_ctx_init(gr_ctx_t ctx, slong nvars, const fmpq_poly_t mu);

/* Coefficient j (0 <= j < n) of an element, an fmpz_mpoly_q over FF_FNX_DATA(ctx)->mctx. */
#define FF_FNX_COEFF(x, j) (((fmpz_mpoly_q_struct*)(x)) + (j))

#endif /* USE_FLINT */
#endif /* FF_FNX_H */
