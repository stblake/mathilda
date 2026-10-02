/*
 * flint_num_bridge.c
 * ------------------
 * FLINT arb/acb-backed numeric evaluation of transcendental special functions.
 * See flint_num_bridge.h. The acb result (a complex ball) is rendered to an
 * Expr through numeric_mpfr_make_complex (the same "Im rounds to 0 -> real
 * leaf" convention the hand-rolled kernels use). Everything is behind
 * USE_FLINT && USE_MPFR (acb -> Expr needs MPFR); the fallback build provides
 * stubs so the module links either way.
 */

#include "flint_num_bridge.h"
#include "symtab.h"
#include "attr.h"
#include "sym_names.h"
#include <stdlib.h>
#include <string.h>

#if defined(USE_FLINT) && defined(USE_MPFR)

#include "numeric_complex.h"   /* numeric_mpfr_make_complex */
#include "numeric.h"           /* numeric_min_inexact_bits  */
#include <gmp.h>
#include <mpfr.h>
#include <flint/fmpz.h>
#include <flint/fmpq.h>
#include <flint/arb.h>
#include <flint/acb.h>
#include <flint/acb_dirichlet.h>
#include <flint/acb_elliptic.h>

/* ------------------------------------------------------------------ */
/*  Expr -> acb  (numeric scalars only)                                */
/* ------------------------------------------------------------------ */

/* A real numeric scalar (Integer / BigInt / Rational / Real / MPFR) -> arb.
 * Set exactly (radius 0) so the ball's midpoint carries the input value; the
 * kernel's own error bounds then govern the result. Returns 0 for anything
 * else (symbolic / Complex / compound), letting the caller bail. */
static int scalar_to_arb(const Expr* e, arb_t out, slong prec) {
    switch (e->type) {
        case EXPR_INTEGER:
            arb_set_si(out, e->data.integer);
            return 1;
        case EXPR_BIGINT: {
            fmpz_t z; fmpz_init(z);
            fmpz_set_mpz(z, e->data.bigint);
            arb_set_fmpz(out, z);
            fmpz_clear(z);
            return 1;
        }
        case EXPR_REAL:
            arb_set_d(out, e->data.real);
            return 1;
        case EXPR_MPFR:
            arf_set_mpfr(arb_midref(out), e->data.mpfr);
            mag_zero(arb_radref(out));
            return 1;
        case EXPR_FUNCTION: {
            const Expr* h = e->data.function.head;
            if (h && h->type == EXPR_SYMBOL
                && strcmp(h->data.symbol.name, "Rational") == 0
                && e->data.function.arg_count == 2) {
                const Expr* a = e->data.function.args[0];
                const Expr* b = e->data.function.args[1];
                if ((a->type != EXPR_INTEGER && a->type != EXPR_BIGINT) ||
                    (b->type != EXPR_INTEGER && b->type != EXPR_BIGINT)) return 0;
                fmpq_t q; fmpq_init(q);
                if (a->type == EXPR_INTEGER) fmpz_set_si(fmpq_numref(q), a->data.integer);
                else                          fmpz_set_mpz(fmpq_numref(q), a->data.bigint);
                if (b->type == EXPR_INTEGER) fmpz_set_si(fmpq_denref(q), b->data.integer);
                else                          fmpz_set_mpz(fmpq_denref(q), b->data.bigint);
                fmpq_canonicalise(q);
                arb_set_fmpq(out, q, prec);
                fmpq_clear(q);
                return 1;
            }
            return 0;
        }
        default:
            return 0;
    }
}

/* A numeric scalar, possibly Complex[re, im], -> acb. */
static int expr_to_acb(const Expr* e, acb_t out, slong prec) {
    if (e->type == EXPR_FUNCTION) {
        const Expr* h = e->data.function.head;
        if (h && h->type == EXPR_SYMBOL
            && strcmp(h->data.symbol.name, "Complex") == 0
            && e->data.function.arg_count == 2) {
            return scalar_to_arb(e->data.function.args[0], acb_realref(out), prec)
                && scalar_to_arb(e->data.function.args[1], acb_imagref(out), prec);
        }
    }
    arb_zero(acb_imagref(out));
    return scalar_to_arb(e, acb_realref(out), prec);
}

/* Extract an fmpz from an integer-like Expr. */
static int expr_to_fmpz(const Expr* e, fmpz_t out) {
    if (e->type == EXPR_INTEGER) { fmpz_set_si(out, e->data.integer); return 1; }
    if (e->type == EXPR_BIGINT)  { fmpz_set_mpz(out, e->data.bigint); return 1; }
    return 0;
}

/* ------------------------------------------------------------------ */
/*  acb -> Expr                                                        */
/* ------------------------------------------------------------------ */

/* Render a finite acb ball to a numeric Expr at `out_bits` precision. NULL if
 * the result is not finite (pole / overflow / indeterminate). */
static Expr* acb_to_expr(const acb_t z, slong out_bits) {
    if (!acb_is_finite(z)) return NULL;
    mpfr_prec_t mp = out_bits > 2 ? (mpfr_prec_t)out_bits : 53;
    mpfr_t re, im;
    mpfr_init2(re, mp);
    mpfr_init2(im, mp);
    arf_get_mpfr(re, arb_midref(acb_realref(z)), MPFR_RNDN);
    arf_get_mpfr(im, arb_midref(acb_imagref(z)), MPFR_RNDN);
    Expr* out = numeric_mpfr_make_complex(re, im);
    mpfr_clear(re);
    mpfr_clear(im);
    return out;
}

/* Target precision from the numeric arguments: the min inexact-bit count over
 * the arguments (an exact arg imposes no constraint), floored at 53. */
static slong pick_out_bits(const Expr* const* args, int n) {
    long best = 0;
    for (int i = 0; i < n; i++) {
        long b = numeric_min_inexact_bits(args[i]);
        if (b > 0 && (best == 0 || b < best)) best = b;
    }
    if (best < 53) best = 53;
    return (slong)best;
}

/* Guard bits added to the working precision above the target. */
#define NB_GUARD 32

/* How many times the working precision may double before we give up asking for
 * `outb` good bits and emit only the bits we actually have. Six doublings is a
 * 64x precision ladder, which covers every conditioning problem seen in
 * practice while keeping the worst case bounded. */
#define NB_MAX_DOUBLINGS 6

/* ------------------------------------------------------------------ */
/*  The evaluation ladder                                              */
/* ------------------------------------------------------------------ */
/* A kernel called once at `outb + NB_GUARD` bits and rendered at `outb` bits is
 * only right when the kernel's own ball stays within NB_GUARD bits. When it does
 * not, rendering the midpoint at `outb` bits prints digits nobody computed. Arb
 * is rigorous and will tell us, so ask.
 *
 * MEASURED, so this is not a speculative guard: N[EllipticPi[-10^30, 1/2], 30]
 * asks for 100 bits and acb_elliptic_pi comes back with 78. Without the ladder
 * the last six of the thirty digits were garbage; with it they are right
 * (1.570796326794897122662117945335*10^-15, mpmath-confirmed).
 *
 * What this is NOT. The ladder cannot see an argument that was ALREADY rounded
 * before it got here, because scalar_to_arb sets every input exactly, radius 0
 * (see its comment) -- so an ill-conditioned function applied to a correctly
 * rounded input gives a TIGHT ball around the wrong value, and no amount of
 * asking the ball will reveal it. That failure mode is real and was the headline
 * bug here (N[EllipticK[99999999999999999/10^17], 20] gave six good digits of
 * twenty; N[Zeta[1 + 1/10^20], 20] gave 2^66 + 1), and it is fixed one layer up,
 * in numeric_plan_working_spec (src/numeric.c), which is what chooses the
 * precision the arguments arrive at. Keep the two straight: this ladder covers
 * the KERNEL's losses, that plan covers the ARGUMENT's.
 *
 * The rules, in order:
 *
 *   finite?            no  -> a genuine pole. Return NULL (stay unevaluated), and
 *                             do NOT retry: doubling cannot make it finite.
 *   rel accuracy       >= outb -> done. This is the first iteration for every
 *                             well-conditioned argument, so the common path pays
 *                             one extra cheap predicate and nothing else.
 *   contains zero?     yes -> relative accuracy is unattainable (the true value
 *                             may BE zero), so retrying is pure cost. Emit the
 *                             midpoint, as this code always did.
 *   otherwise          climb; on exhausting the ladder emit only the bits the
 *                             ball actually has, which is fewer digits and all
 *                             of them correct.
 */

/* Significant bits worth rendering from `z` when `outb` were asked for, or -1 to
 * mean "not finite -- bail". */
static slong nb_emit_bits(const acb_t z, slong outb) {
    if (!acb_is_finite(z)) return -1;
    slong acc = acb_rel_accuracy_bits(z);
    if (acc >= outb || acb_contains_zero(z)) return outb;
    return acc > 2 ? acc : outb;   /* acc <= 2: nothing to report; keep old shape */
}

/* True when another doubling could plausibly help. */
static int nb_should_retry(const acb_t z, slong outb, int attempt) {
    if (attempt >= NB_MAX_DOUBLINGS) return 0;
    if (!acb_is_finite(z)) return 0;              /* pole: hopeless */
    if (acb_contains_zero(z)) return 0;           /* relative target unattainable */
    return acb_rel_accuracy_bits(z) < outb;
}

/* Arb kernel shapes, by arity. The incomplete elliptic kernels carry an extra
 * `times_pi` flag; static adapters below absorb it so they fit these types. */
typedef void (*nb_fn1)(acb_t, const acb_t, slong);
typedef void (*nb_fn2)(acb_t, const acb_t, const acb_t, slong);
typedef void (*nb_fn3)(acb_t, const acb_t, const acb_t, const acb_t, slong);

/* Evaluate `f` at `args`, climbing the precision ladder. The arguments are
 * re-converted on every attempt because scalar_to_arb ROUNDS a Rational at
 * `prec` -- a retry that reused the first conversion would raise the working
 * precision around an input that is still only accurate to the old one. */
static Expr* nb_eval(const Expr* const* args, int n, slong outb,
                     nb_fn1 f1, nb_fn2 f2, nb_fn3 f3) {
    acb_t A[3], R;
    for (int i = 0; i < n; i++) acb_init(A[i]);
    acb_init(R);
    Expr* out = NULL;
    for (int attempt = 0; ; attempt++) {
        slong wp = (outb + NB_GUARD) << attempt;
        int ok = 1;
        for (int i = 0; i < n && ok; i++) ok = expr_to_acb(args[i], A[i], wp);
        if (!ok) break;                            /* non-numeric: decline */
        if (n == 1)      f1(R, A[0], wp);
        else if (n == 2) f2(R, A[0], A[1], wp);
        else             f3(R, A[0], A[1], A[2], wp);
        if (nb_should_retry(R, outb, attempt)) continue;
        slong bits = nb_emit_bits(R, outb);
        if (bits > 0) out = acb_to_expr(R, bits);
        break;
    }
    for (int i = 0; i < n; i++) acb_clear(A[i]);
    acb_clear(R);
    return out;
}

static Expr* nb_eval1(const Expr* a, nb_fn1 f) {
    const Expr* args[1] = { a };
    return nb_eval(args, 1, pick_out_bits(args, 1), f, NULL, NULL);
}

static Expr* nb_eval2(const Expr* a, const Expr* b, nb_fn2 f) {
    const Expr* args[2] = { a, b };
    return nb_eval(args, 2, pick_out_bits(args, 2), NULL, f, NULL);
}

static Expr* nb_eval3(const Expr* a, const Expr* b, const Expr* c, nb_fn3 f) {
    const Expr* args[3] = { a, b, c };
    return nb_eval(args, 3, pick_out_bits(args, 3), NULL, NULL, f);
}

/* ------------------------------------------------------------------ */
/*  Kernels                                                            */
/* ------------------------------------------------------------------ */

Expr* flint_num_zeta(const Expr* s) {
    return nb_eval1(s, acb_dirichlet_zeta);
}

Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a) {
    return nb_eval2(s, a, acb_dirichlet_hurwitz);
}

Expr* flint_num_polygamma(const Expr* n, const Expr* z) {
    return nb_eval2(n, z, acb_polygamma);
}

/* StieltjesGamma[n] or StieltjesGamma[n, a]; n a non-negative integer. The fmpz
 * order does not fit the nb_fn* shapes, so this one climbs the ladder by hand.
 * `a` governs the precision; the default a = 1 is exact, hence the 53-bit floor. */
Expr* flint_num_stieltjes(const Expr* n, const Expr* a) {
    fmpz_t N; fmpz_init(N);
    if (!expr_to_fmpz(n, N) || fmpz_sgn(N) < 0) { fmpz_clear(N); return NULL; }
    const Expr* args[1] = { a ? a : n };
    slong outb = a ? pick_out_bits(args, 1) : 53;
    acb_t A, R; acb_init(A); acb_init(R);
    Expr* out = NULL;
    for (int attempt = 0; ; attempt++) {
        slong wp = (outb + NB_GUARD) << attempt;
        if (a && !expr_to_acb(a, A, wp)) break;
        if (!a) acb_one(A);
        acb_dirichlet_stieltjes(R, N, A, wp);
        if (nb_should_retry(R, outb, attempt)) continue;
        slong bits = nb_emit_bits(R, outb);
        if (bits > 0) out = acb_to_expr(R, bits);
        break;
    }
    acb_clear(A); acb_clear(R);
    fmpz_clear(N);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Legendre elliptic integrals                                        */
/* ------------------------------------------------------------------ */
/* Arb's acb_elliptic_* already use the PARAMETER convention m = k^2 and
 * Mathematica's branch placement, so these are straight pass-throughs. The
 * `times_pi` flag is 0 throughout: phi arrives in radians, not as a multiple of
 * Pi, and the adapters below pin it so the kernels fit the nb_fn* shapes. */

static void nb_elliptic_f(acb_t r, const acb_t phi, const acb_t m, slong prec) {
    acb_elliptic_f(r, phi, m, 0, prec);
}

static void nb_elliptic_e_inc(acb_t r, const acb_t phi, const acb_t m, slong prec) {
    acb_elliptic_e_inc(r, phi, m, 0, prec);
}

static void nb_elliptic_pi_inc(acb_t r, const acb_t n, const acb_t phi,
                               const acb_t m, slong prec) {
    acb_elliptic_pi_inc(r, n, phi, m, 0, prec);
}

Expr* flint_num_elliptic_k(const Expr* m) {
    return nb_eval1(m, acb_elliptic_k);
}

Expr* flint_num_elliptic_e(const Expr* m) {
    return nb_eval1(m, acb_elliptic_e);
}

Expr* flint_num_elliptic_f(const Expr* phi, const Expr* m) {
    return nb_eval2(phi, m, nb_elliptic_f);
}

Expr* flint_num_elliptic_e_inc(const Expr* phi, const Expr* m) {
    return nb_eval2(phi, m, nb_elliptic_e_inc);
}

Expr* flint_num_elliptic_pi(const Expr* n, const Expr* m) {
    return nb_eval2(n, m, acb_elliptic_pi);
}

Expr* flint_num_elliptic_pi_inc(const Expr* n, const Expr* phi, const Expr* m) {
    return nb_eval3(n, phi, m, nb_elliptic_pi_inc);
}

/* ------------------------------------------------------------------ */
/*  FLINT` context builtins                                            */
/* ------------------------------------------------------------------ */

static Expr* builtin_flint_zeta(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    return flint_num_zeta(res->data.function.args[0]);
}

static Expr* builtin_flint_hurwitzzeta(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return flint_num_hurwitz_zeta(res->data.function.args[0],
                                  res->data.function.args[1]);
}

static Expr* builtin_flint_polygamma(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    return flint_num_polygamma(res->data.function.args[0],
                               res->data.function.args[1]);
}

static Expr* builtin_flint_stieltjes(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t n = res->data.function.arg_count;
    if (n == 1) return flint_num_stieltjes(res->data.function.args[0], NULL);
    if (n == 2) return flint_num_stieltjes(res->data.function.args[0],
                                           res->data.function.args[1]);
    return NULL;
}

void flint_num_bridge_init(void) {
    symtab_add_builtin(SYM_FLINT_Zeta, builtin_flint_zeta);
    symtab_get_def(SYM_FLINT_Zeta)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_Zeta,
        "FLINT`Zeta[s] gives the numeric value of the Riemann zeta function at "
        "the numeric argument s (real or complex), computed to the precision "
        "of s (machine precision for exact s) via FLINT's rigorous acb "
        "arithmetic (acb_dirichlet_zeta). Unevaluated for symbolic s or at the "
        "pole s = 1.");

    symtab_add_builtin(SYM_FLINT_HurwitzZeta, builtin_flint_hurwitzzeta);
    symtab_get_def(SYM_FLINT_HurwitzZeta)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_HurwitzZeta,
        "FLINT`HurwitzZeta[s, a] gives the numeric value of the Hurwitz zeta "
        "function via FLINT (acb_dirichlet_hurwitz), to the precision of the "
        "arguments. Unevaluated for symbolic arguments or at a pole.");

    symtab_add_builtin(SYM_FLINT_PolyGamma, builtin_flint_polygamma);
    symtab_get_def(SYM_FLINT_PolyGamma)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_PolyGamma,
        "FLINT`PolyGamma[n, z] gives the numeric value of the n-th derivative "
        "of the digamma function (n = 0 is digamma) via FLINT "
        "(acb_polygamma), to the precision of the arguments. Unevaluated for "
        "symbolic arguments or at a pole.");

    symtab_add_builtin(SYM_FLINT_StieltjesGamma, builtin_flint_stieltjes);
    symtab_get_def(SYM_FLINT_StieltjesGamma)->attributes |= ATTR_PROTECTED;
    symtab_set_docstring(SYM_FLINT_StieltjesGamma,
        "FLINT`StieltjesGamma[n] and FLINT`StieltjesGamma[n, a] give the "
        "numeric value of the n-th Stieltjes constant (generalized, at a) for "
        "a non-negative integer n, via FLINT (acb_dirichlet_stieltjes). "
        "Unevaluated for negative or non-integer n.");
}

#else /* !(USE_FLINT && USE_MPFR) */

Expr* flint_num_zeta(const Expr* s) { (void)s; return NULL; }
Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a) { (void)s; (void)a; return NULL; }
Expr* flint_num_polygamma(const Expr* n, const Expr* z) { (void)n; (void)z; return NULL; }
Expr* flint_num_stieltjes(const Expr* n, const Expr* a) { (void)n; (void)a; return NULL; }
Expr* flint_num_elliptic_k(const Expr* m) { (void)m; return NULL; }
Expr* flint_num_elliptic_e(const Expr* m) { (void)m; return NULL; }
Expr* flint_num_elliptic_f(const Expr* phi, const Expr* m) { (void)phi; (void)m; return NULL; }
Expr* flint_num_elliptic_e_inc(const Expr* phi, const Expr* m) { (void)phi; (void)m; return NULL; }
Expr* flint_num_elliptic_pi(const Expr* n, const Expr* m) { (void)n; (void)m; return NULL; }
Expr* flint_num_elliptic_pi_inc(const Expr* n, const Expr* phi, const Expr* m) {
    (void)n; (void)phi; (void)m; return NULL;
}
void  flint_num_bridge_init(void) { /* no FLINT/MPFR: nothing to register */ }

#endif
