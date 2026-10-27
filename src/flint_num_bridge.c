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
#include "numeric.h"           /* arg_is_inexact, numeric_lower_mpfr_to_machine */

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
#include <flint/dirichlet.h>
#include <flint/ulong_extras.h>
#include <flint/acb_elliptic.h>
#include <flint/acb_hypgeom.h>
#include <flint/arb_hypgeom.h>
#include <flint/acb_poly.h>

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

/* True when the ball's absolute radius (either part) is at least 2^-outb, i.e.
 * the ball is not tight around zero in absolute terms either. */
static int nb_radius_not_tiny(const acb_t z, slong outb) {
    return mag_cmp_2exp_si(arb_radref(acb_realref(z)), -outb) >= 0 ||
           mag_cmp_2exp_si(arb_radref(acb_imagref(z)), -outb) >= 0;
}

/* Significant bits worth rendering from `z` when `outb` were asked for, or -1 to
 * mean "not finite -- bail". */
static slong nb_emit_bits(const acb_t z, slong outb) {
    if (!acb_is_finite(z)) return -1;
    slong acc = acb_rel_accuracy_bits(z);
    if (acc >= outb) return outb;
    if (acb_contains_zero(z))   /* tight: emit the midpoint; wide: no digits known */
        return nb_radius_not_tiny(z, outb) ? -1 : outb;
    return acc > 2 ? acc : outb;   /* acc <= 2: nothing to report; keep old shape */
}

/* True when another doubling could plausibly help. */
static int nb_should_retry(const acb_t z, slong outb, int attempt) {
    if (attempt >= NB_MAX_DOUBLINGS) return 0;
    if (!acb_is_finite(z)) return 0;              /* pole: hopeless */
    /* A ball containing zero: relative accuracy is unattainable when the true
     * value may BE zero, so a TIGHT ball (radius < 2^-outb) is accepted.  A
     * WIDE one is a cancellation loss, not a zero -- HurwitzZeta[-38.5`25,
     * 1/2`25] came back from acb_dirichlet_hurwitz at 116 bits as
     * 2.5e24 +/- 1.1e25 (true value -1.36e14), and emitting that midpoint is
     * how LerchPhi's |z| > 1 continuation went wrong at 25 digits.  Climb. */
    if (acb_contains_zero(z)) return nb_radius_not_tiny(z, outb);
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

/* PolyGamma[nu, z] for a NON-INTEGER order, in the reference's convention:
 * the Riemann-Liouville derivative of order nu, from 0, of psi(z) + 1/z (which
 * is regular there), with the pole term 1/z continued analytically in nu:
 *
 *   psi^(nu)(z) = z^-nu/Gamma(1-nu) [ -gamma + Sum_{k>=1} (1 - 2F1(1,1;1-nu;-z/k))/k ]
 *                 - z^(-1-nu)/Gamma(-nu) (Log z - psi(-nu) - gamma)
 *
 * using RL[(k+z)^-1] = z^-nu/(k Gamma(1-nu)) 2F1(1,1;1-nu;-z/k). The sum runs
 * k = 1..K directly and the tail k > K through its Taylor expansion in z,
 *
 *   Sum_{k>K} (1/k - 1/(k+z)) -> -Sum_{n>=1} n!/(1-nu)_n (-z)^n zeta(n+1, K+1),
 *
 * which converges geometrically (ratio |z|/(K+1) <= 1/2 for K = 2|z| + 2). It
 * reduces to the ordinary polygamma as nu approaches an integer, and agreed with
 * the reference to the 20 digits it printed at nu = -1/2, 1/2, 3/2 and z = 1, 2,
 * -3/2, 2 + I. |z| > 5000 declines (the direct sum would cost 10^4 2F1s). */
static void nb_polygamma_frac(acb_t res, const acb_t nu, const acb_t z_in, slong prec) {
    /* A local copy: GCC's -Wstringop-overread misreads arb's inline helpers
     * applied straight to the const parameter. */
    acb_t z; acb_init(z); acb_set(z, z_in);
    double zr = arf_get_d(arb_midref(acb_realref(z)), ARF_RND_NEAR);
    double zi = arf_get_d(arb_midref(acb_imagref(z)), ARF_RND_NEAR);
    double az = sqrt(zr * zr + zi * zi);
    if (!(az <= 5000.0)) { acb_indeterminate(res); acb_clear(z); return; }
    slong K = (slong)(2.0 * az) + 2;
    double ratio = az / (double)(K + 1);
    slong nmax = (slong)((double)(prec + 20) / -log2(ratio > 1e-300 ? ratio : 1e-300)) + 20;
    if (nmax > 4 * prec + 100) nmax = 4 * prec + 100;

    acb_t one, t, g1, zn, c, eul, S, sum, w, term, coef, pw, mz, hz, s1, hk, a;
    acb_init(one); acb_init(t); acb_init(g1); acb_init(zn); acb_init(c);
    acb_init(eul); acb_init(S); acb_init(sum); acb_init(w); acb_init(term);
    acb_init(coef); acb_init(pw); acb_init(mz); acb_init(hz); acb_init(s1);
    acb_init(hk); acb_init(a);

    acb_one(one);
    acb_sub(t, one, nu, prec);                  /* t = 1 - nu */
    acb_rgamma(g1, t, prec);                    /* 1/Gamma(1 - nu) */
    acb_neg(a, nu);
    acb_pow(zn, z, a, prec);                    /* z^-nu */
    acb_mul(c, zn, g1, prec);                   /* c = z^-nu/Gamma(1-nu) */
    arb_const_euler(acb_realref(eul), prec);

    /* -gamma + the direct terms k = 1..K. */
    acb_neg(sum, eul);
    for (slong k = 1; k <= K; k++) {
        acb_div_si(w, z, -k, prec);
        acb_hypgeom_2f1(term, one, one, t, w, 0, prec);
        acb_sub(term, one, term, prec);
        acb_div_si(term, term, k, prec);
        acb_add(sum, sum, term, prec);
    }
    /* The tail k > K. */
    acb_neg(mz, z);
    acb_one(coef); acb_one(pw); acb_zero(s1);
    acb_set_si(hk, K + 1);
    for (slong n = 1; n <= nmax; n++) {
        acb_add_si(w, t, n - 1, prec);          /* (1-nu)_n / (1-nu)_{n-1} */
        acb_mul_si(coef, coef, n, prec);
        acb_div(coef, coef, w, prec);
        acb_mul(pw, pw, mz, prec);
        acb_set_si(w, n + 1);
        acb_dirichlet_hurwitz(hz, w, hk, prec);
        acb_mul(term, coef, pw, prec);
        acb_mul(term, term, hz, prec);
        acb_add(s1, s1, term, prec);
    }
    acb_sub(sum, sum, s1, prec);
    acb_mul(S, c, sum, prec);

    /* The continued pole term. */
    acb_rgamma(g1, a, prec);                    /* 1/Gamma(-nu) */
    acb_digamma(w, a, prec);
    acb_log(term, z, prec);
    acb_sub(term, term, w, prec);
    acb_sub(term, term, eul, prec);
    acb_sub_si(w, a, 1, prec);
    acb_pow(hz, z, w, prec);                    /* z^(-1-nu) */
    acb_mul(term, term, hz, prec);
    acb_mul(term, term, g1, prec);
    acb_sub(res, S, term, prec);

    acb_clear(one); acb_clear(t); acb_clear(g1); acb_clear(zn); acb_clear(c);
    acb_clear(eul); acb_clear(S); acb_clear(sum); acb_clear(w); acb_clear(term);
    acb_clear(coef); acb_clear(pw); acb_clear(mz); acb_clear(hz); acb_clear(s1);
    acb_clear(hk); acb_clear(a); acb_clear(z);
}

Expr* flint_num_polygamma_frac(const Expr* nu, const Expr* z) {
    return nb_eval2(nu, z, nb_polygamma_frac);
}

/* (a; q)_inf = prod_{k>=0} (1 - a q^k) for |q| < 1, at `prec`: the factors are
 * multiplied until |a q^k| < 2^-(prec + 10), which bounds the rest of the
 * product to within that of 1. Indeterminate (so the caller declines) when
 * |q| >= 1 or 10^6 factors do not reach it. */
static void nb_qpoch_inf(acb_t res, const acb_t a, const acb_t q, slong prec) {
    acb_t aqk, f, qq;
    acb_init(aqk); acb_init(f); acb_init(qq);
    acb_set(qq, q);
    mag_t mq; mag_init(mq);
    acb_get_mag(mq, qq);
    if (mag_cmp_2exp_si(mq, 0) >= 0) { acb_indeterminate(res); goto done; }
    acb_one(res);
    acb_set(aqk, a);
    for (slong k = 0; ; k++) {
        if (k >= 1000000) { acb_indeterminate(res); break; }
        acb_sub_si(f, aqk, 1, prec);
        acb_neg(f, f);                       /* 1 - a q^k */
        acb_mul(res, res, f, prec);
        acb_mul(aqk, aqk, qq, prec);
        mag_t m; mag_init(m);
        acb_get_mag(m, aqk);
        int small = mag_cmp_2exp_si(m, -(prec + 10)) < 0;
        mag_clear(m);
        if (small) break;
    }
done:
    mag_clear(mq);
    acb_clear(aqk); acb_clear(f); acb_clear(qq);
}

/* (a; q)_n = (a; q)_inf / (a q^n; q)_inf, for any complex n (|q| < 1). */
static void nb_qpoch_n(acb_t res, const acb_t a, const acb_t q, const acb_t n, slong prec) {
    acb_t qn, den;
    acb_init(qn); acb_init(den);
    acb_pow(qn, q, n, prec);
    acb_mul(qn, qn, a, prec);
    nb_qpoch_inf(den, qn, q, prec);
    nb_qpoch_inf(res, a, q, prec);
    acb_div(res, res, den, prec);
    acb_clear(qn); acb_clear(den);
}

Expr* flint_num_qpochhammer(const Expr* a, const Expr* q, const Expr* n) {
    return n ? nb_eval3(a, q, n, nb_qpoch_n) : nb_eval2(a, q, nb_qpoch_inf);
}

/* The zero of J_nu nearest `seed` (a machine zero), refined by Newton in arb at
 * `bits` (+ guard) bits: x <- x - J_nu(x) / J_nu'(x), J' = (J_{nu-1} - J_{nu+1})/2,
 * with the working precision doubling from 64 bits so each step starts from an
 * iterate already correct to the previous precision. */
Expr* flint_num_bessel_j_zero(const Expr* nu, double seed, long bits) {
    slong target = bits + NB_GUARD;
    acb_t N, Nm, Np, X, F, A, B;
    acb_init(N); acb_init(Nm); acb_init(Np); acb_init(X); acb_init(F);
    acb_init(A); acb_init(B);
    Expr* out = NULL;
    if (!expr_to_acb(nu, N, target)) goto done;
    acb_set_d(X, seed);
    for (slong wp = 64; ; wp *= 2) {
        if (wp > target) wp = target;
        acb_sub_si(Nm, N, 1, wp);
        acb_add_si(Np, N, 1, wp);
        for (int it = 0; it < 4; it++) {
            acb_hypgeom_bessel_j(F, N, X, wp);
            acb_hypgeom_bessel_j(A, Nm, X, wp);
            acb_hypgeom_bessel_j(B, Np, X, wp);
            acb_sub(A, A, B, wp);
            acb_mul_2exp_si(A, A, -1);
            acb_div(F, F, A, wp);
            acb_sub(X, X, F, wp);
            acb_get_mid(X, X);          /* Newton: carry the midpoint */
        }
        if (wp == target) break;
    }
    if (acb_is_finite(X)) out = acb_to_expr(X, bits);
done:
    acb_clear(N); acb_clear(Nm); acb_clear(Np); acb_clear(X); acb_clear(F);
    acb_clear(A); acb_clear(B);
    return out;
}

/* Nielsen's generalized polylogarithm S_{n,p}(z), a positive integer p and
 * |z| < 1, from (-Log(1 - w))^p / p! = Sum_k c(k, p) w^k / k! (c the unsigned
 * Stirling numbers of the first kind) integrated against the Mellin form
 * S_{n,p}(z) = (-1)^p / (Gamma(n) p!) Int_0^inf u^(n-1) Log^p(1 - z e^-u) du:
 *
 *   S_{n,p}(z) = Sum_{k>=p} (c(k, p) / k!) z^k / k^n.
 *
 * e(k, j) = c(k, j)/k! is carried for j = 0..p by
 * e(k+1, j) = (k e(k, j) + e(k, j-1)) / (k+1). The terms decay like
 * |z|^k (Log k)^(p-1) / k^(n+1); the sum stops once |term| < 2^-(prec+10) past
 * the peak, and declines (Indeterminate) beyond 10^6 terms or for |z| >= 1,
 * where the series does not converge. It is the small-|z| method of
 * nb_nielsen_gen for an integer p. */
static void nb_nielsen(acb_t res, const acb_t n, slong p, const acb_t z_in, slong prec) {
    acb_t z; acb_init(z); acb_set(z, z_in);
    mag_t mz; mag_init(mz); acb_get_mag(mz, z);
    double az = mag_get_d(mz);
    mag_clear(mz);
    if (!(az < 1.0) || p < 1 || p > 1000) { acb_indeterminate(res); acb_clear(z); return; }
    slong nmax = az == 0.0 ? p + 1
               : (slong)((double)(prec + 20) * 0.6931471805599453 / -log(az)) + 4 * p + 100;
    if (nmax > 1000000) { acb_indeterminate(res); acb_clear(z); return; }

    arb_struct* e = (arb_struct*)flint_malloc(sizeof(arb_struct) * (size_t)(p + 1));
    for (slong j = 0; j <= p; j++) arb_init(e + j);
    acb_t zk, term, kn, kk;
    acb_init(zk); acb_init(term); acb_init(kn); acb_init(kk);
    arb_t tmp; arb_init(tmp);
    acb_zero(res);
    arb_one(e + 0);                         /* k = 0: e(0,0) = 1 */
    acb_one(zk);
    for (slong k = 0; k < nmax; k++) {
        /* advance to k + 1 */
        for (slong j = p; j >= 1; j--) {
            arb_mul_si(tmp, e + j, k, prec);
            arb_add(tmp, tmp, e + j - 1, prec);
            arb_div_si(e + j, tmp, k + 1, prec);
        }
        arb_mul_si(e + 0, e + 0, k, prec);  /* e(k+1, 0) = 0 for k >= 0 */
        arb_div_si(e + 0, e + 0, k + 1, prec);
        acb_mul(zk, zk, z, prec);
        slong kp = k + 1;
        if (kp < p) continue;
        acb_set_si(kk, kp);
        acb_pow(kn, kk, n, prec);
        acb_div(term, zk, kn, prec);
        acb_mul_arb(term, term, e + p, prec);
        acb_add(res, res, term, prec);
        if (kp > 4 * p + 20) {
            mag_t mt; mag_init(mt); acb_get_mag(mt, term);
            int tiny = mag_cmp_2exp_si(mt, -(prec + 10)) < 0;
            mag_clear(mt);
            if (tiny) break;
        }
    }
    for (slong j = 0; j <= p; j++) arb_clear(e + j);
    flint_free(e);
    acb_clear(zk); acb_clear(term); acb_clear(kn); acb_clear(kk);
    arb_clear(tmp);
    acb_clear(z);
}

/* S_{n,p}(z) for a non-integer (or complex) p and |z| < 1. With g(w) =
 * -Log(1 - w)/w = Sum_j w^j/(j+1), (-Log(1 - w))^p = w^p g(w)^p, and the
 * coefficients c_k of g^p follow J. C. P. Miller's power recurrence
 * c_k = (1/k) Sum_{j=1}^{k} ((p + 1) j - k) g_j c_{k-j}, so term by term in the
 * Mellin form
 *
 *   S_{n,p}(z) = (1/Gamma(p + 1)) Sum_{k>=0} c_k z^(p+k) / (p + k)^n.
 *
 * O(K^2) in the K ~ prec / -log2|z| terms, so it is the method for small |z|
 * or Re n <= 0 (where the quadrature below does not converge); K is capped at
 * 4000. */
static void nb_nielsen_series_gen(acb_t res, const acb_t n, const acb_t p,
                                  const acb_t z, slong prec) {
    mag_t mz; mag_init(mz); acb_get_mag(mz, z);
    double az = mag_get_d(mz);
    mag_clear(mz);
    if (!(az < 1.0)) { acb_indeterminate(res); return; }
    slong K = az == 0.0 ? 1 : (slong)((double)(prec + 20) / -log2(az)) + 10;
    if (K > 4000) { acb_indeterminate(res); return; }
    acb_struct* c = _acb_vec_init(K + 1);
    acb_t s, t, zp, zk, pk, w, p1;
    acb_init(s); acb_init(t); acb_init(zp); acb_init(zk); acb_init(pk);
    acb_init(w); acb_init(p1);
    acb_add_ui(p1, p, 1, prec);                     /* p + 1 */
    acb_one(c + 0);
    for (slong k = 1; k <= K; k++) {
        acb_zero(s);
        for (slong j = 1; j <= k; j++) {
            acb_mul_si(t, p1, j, prec);                 /* (p + 1) j - k */
            acb_sub_si(t, t, k, prec);
            acb_div_si(t, t, j + 1, prec);              /* g_j = 1/(j + 1) */
            acb_addmul(s, t, c + k - j, prec);
        }
        acb_div_si(c + k, s, k, prec);
    }
    acb_pow(zp, z, p, prec);                        /* z^p, principal */
    acb_one(zk);
    acb_zero(res);
    for (slong k = 0; k <= K; k++) {
        acb_add_si(pk, p, k, prec);
        acb_pow(w, pk, n, prec);                    /* (p + k)^n */
        acb_mul(t, c + k, zk, prec);
        acb_div(t, t, w, prec);
        acb_add(res, res, t, prec);
        acb_mul(zk, zk, z, prec);
    }
    acb_mul(res, res, zp, prec);
    acb_rgamma(t, p1, prec);
    acb_mul(res, res, t, prec);
    _acb_vec_clear(c, K + 1);
    acb_clear(s); acb_clear(t); acb_clear(zp); acb_clear(zk); acb_clear(pk);
    acb_clear(w); acb_clear(p1);
}

/* Tanh-sinh (double-exponential) quadrature of f over [0, 1], for integrands
 * with integrable singularities at either end. f receives the node t and its
 * complement 1 - t, each computed directly (1/(1 + e^-2s) and e^-2s/(1 + e^-2s))
 * so a logarithm of either keeps full relative accuracy at the endpoints. The
 * step halves until two successive sums agree to 2^-(prec+4); the difference
 * of the last two is added to the result's radius, so the caller's precision
 * ladder sees an unconverged sum as an inaccurate one. Indeterminate when no
 * level converges. */
typedef void (*nb_qfn)(acb_t out, const arb_t t, const arb_t omt, const void* ctx, slong prec);

static void nb_tanh_sinh(acb_t res, nb_qfn f, const void* ctx, slong prec) {
    slong wp = prec + 40;   /* the node sum carries ~2^16 terms' rounding */
    double khmax = asinh((double)(wp + 60) * 0.2206356) + 0.5;  /* ln2/pi */
    int maxlevel = 6 + (int)(log2((double)prec));
    acb_t S, prev, v, d;
    arb_t pi, h, kh, sh, ch, s2, e, onepe, t, omt, wgt;
    acb_init(S); acb_init(prev); acb_init(v); acb_init(d);
    arb_init(pi); arb_init(h); arb_init(kh); arb_init(sh); arb_init(ch);
    arb_init(s2); arb_init(e); arb_init(onepe); arb_init(t); arb_init(omt); arb_init(wgt);
    arb_const_pi(pi, wp);
    int done = 0;
    for (int level = 2; level <= maxlevel && !done; level++) {
        arb_one(h); arb_mul_2exp_si(h, h, -level);
        /* the centre node t = 1/2 */
        arb_set_d(t, 0.5); arb_set_d(omt, 0.5);
        f(v, t, omt, ctx, wp);
        arb_mul(wgt, h, pi, wp); arb_mul_2exp_si(wgt, wgt, -2);    /* h pi/4 */
        acb_mul_arb(S, v, wgt, wp);
        slong kmax = (slong)(khmax * (double)(1 << level)) + 1;
        for (slong k = 1; k <= kmax; k++) {
            arb_mul_si(kh, h, k, wp);
            arb_sinh_cosh(sh, ch, kh, wp);
            arb_mul(s2, sh, pi, wp);                               /* 2s = pi sinh(kh) */
            arb_neg(e, s2); arb_exp(e, e, wp);                     /* e^-2s */
            arb_add_ui(onepe, e, 1, wp);
            arb_inv(t, onepe, wp);                                 /* 1/(1 + e^-2s) */
            arb_mul(omt, e, t, wp);                                /* e^-2s/(1 + e^-2s) */
            /* w = h pi cosh(kh) e^-2s / (1 + e^-2s)^2 */
            arb_mul(wgt, h, pi, wp); arb_mul(wgt, wgt, ch, wp);
            arb_mul(wgt, wgt, omt, wp); arb_mul(wgt, wgt, t, wp);
            f(v, t, omt, ctx, wp);                                 /* the node near 1 */
            if (acb_is_finite(v)) acb_addmul_arb(S, v, wgt, wp);
            f(v, omt, t, ctx, wp);                                 /* its mirror near 0 */
            if (acb_is_finite(v)) acb_addmul_arb(S, v, wgt, wp);
            if (arb_is_zero(omt)) break;
        }
        if (level > 2) {
            /* Midpoints only: the balls' accumulated rounding radius is not the
             * discretisation error, and would keep this test from ever passing. */
            acb_get_mid(d, S);
            acb_sub(d, d, prev, wp);
            acb_get_mid(d, d);
            mag_t md, ms; mag_init(md); mag_init(ms);
            acb_get_mag(md, d); acb_get_mag(ms, S);
            mag_mul_2exp_si(ms, ms, -(prec + 4));
            if (mag_cmp(md, ms) <= 0 && acb_is_finite(S)) {
                acb_get_mid(res, S);
                acb_add_error_mag(res, md);
                done = 1;
            }
            mag_clear(md); mag_clear(ms);
        }
        acb_get_mid(prev, S);
    }
    if (!done) acb_indeterminate(res);
    acb_clear(S); acb_clear(prev); acb_clear(v); acb_clear(d);
    arb_clear(pi); arb_clear(h); arb_clear(kh); arb_clear(sh); arb_clear(ch);
    arb_clear(s2); arb_clear(e); arb_clear(onepe); arb_clear(t); arb_clear(omt); arb_clear(wgt);
}

/* The integrand of S_{n,p}(z) = 1/(Gamma(n) Gamma(p+1)) Int_0^1
 * (-Log t)^(n-1) (-Log(1 - z t))^p dt/t on one of three pieces of [0, 1]:
 * the whole interval, or, for a real z > 1 whose branch point t = 1/z lies on
 * it, [0, 1/z] and [1/z, 1]. Each piece is mapped to tau in [0, 1] so that the
 * quantity that vanishes at its ends is formed without cancellation:
 *   whole:  t = tau,                 1 - z t = (1 - tau) - (z - 1) tau
 *   left:   t = c tau (c = 1/z),     1 - z t = 1 - tau
 *   right:  t = c + (1 - c) tau,     1 - z t = (1 - z) tau   (negative: Log
 *           takes the principal + I Pi, the reference's side of the cut). */
typedef struct {
    acb_t n1, p, z, zm1;   /* n - 1, p, z, z - 1 */
    arb_t c, omc;          /* the split point 1/z and 1 - 1/z */
    int piece;             /* 0 whole, 1 left, 2 right */
    int p_int;             /* p an exact small integer: acb_pow_si */
    slong p_si;
} nb_nielsen_ctx;

static void nb_nielsen_integrand(acb_t out, const arb_t tau, const arb_t omtau,
                                 const void* vctx, slong prec) {
    const nb_nielsen_ctx* C = (const nb_nielsen_ctx*)vctx;
    arb_t t, omt, r;
    acb_t w, lw, lt;
    arb_init(t); arb_init(omt); arb_init(r);
    acb_init(w); acb_init(lw); acb_init(lt);
    if (C->piece == 0) {
        arb_set(t, tau); arb_set(omt, omtau);
        acb_mul_arb(w, C->zm1, tau, prec);
        acb_neg(w, w);
        acb_add_arb(w, w, omtau, prec);
    } else if (C->piece == 1) {
        arb_mul(t, C->c, tau, prec);
        arb_mul(omt, C->c, tau, prec); arb_sub_ui(omt, omt, 1, prec); arb_neg(omt, omt);
        acb_set_arb(w, omtau);
    } else {
        arb_mul(t, C->omc, tau, prec); arb_add(t, t, C->c, prec);
        arb_mul(omt, C->omc, omtau, prec);
        acb_mul_arb(w, C->zm1, tau, prec);
        acb_neg(w, w);
    }
    /* -Log t = Log1p((1 - t)/t), accurate as t -> 1 */
    arb_div(r, omt, t, prec);
    arb_log1p(r, r, prec);
    acb_set_arb(lt, r);
    /* -Log(1 - u) for the u = z t (or tau) that vanishes at the left end,
     * through log1p: Log of 1 - tiny cancels, and a fractional power of the
     * resulting ball about 0 is worthless. */
    if (C->piece != 2) {
        acb_t u; acb_init(u);
        if (C->piece == 0) acb_mul_arb(u, C->z, t, prec);
        else               acb_set_arb(u, tau);
        mag_t mu; mag_init(mu); acb_get_mag(mu, u);
        if (mag_cmp_2exp_si(mu, -1) < 0) {
            acb_neg(u, u);
            acb_log1p(lw, u, prec);
        } else {
            acb_log(lw, w, prec);
        }
        mag_clear(mu); acb_clear(u);
    } else {
        acb_log(lw, w, prec);
    }
    acb_neg(lw, lw);
    if (C->p_int) acb_pow_si(lw, lw, C->p_si, prec);
    else          acb_pow(lw, lw, C->p, prec);
    acb_pow(lt, lt, C->n1, prec);
    acb_mul(out, lt, lw, prec);
    acb_div_arb(out, out, t, prec);
    if (C->piece == 1) acb_mul_arb(out, out, C->c, prec);
    if (C->piece == 2) acb_mul_arb(out, out, C->omc, prec);
    arb_clear(t); arb_clear(omt); arb_clear(r);
    acb_clear(w); acb_clear(lw); acb_clear(lt);
}

/* S_{n,p}(z) by quadrature of its Mellin-type integral; Re n > 0. */
static void nb_nielsen_quad(acb_t res, const acb_t n, const acb_t p, const acb_t z, slong prec) {
    nb_nielsen_ctx C;
    acb_init(C.n1); acb_init(C.p); acb_init(C.z); acb_init(C.zm1);
    arb_init(C.c); arb_init(C.omc);
    acb_sub_ui(C.n1, n, 1, prec);
    acb_set(C.p, p);
    acb_set(C.z, z);
    acb_sub_ui(C.zm1, z, 1, prec);
    C.p_int = acb_is_int(p) && arb_is_exact(acb_realref(p)) &&
              arf_cmpabs_2exp_si(arb_midref(acb_realref(p)), 20) < 0;
    C.p_si = C.p_int ? arf_get_si(arb_midref(acb_realref(p)), ARF_RND_NEAR) : 0;
    acb_t part, g;
    acb_init(part); acb_init(g);
    int split = acb_is_real(z) && arf_cmp_si(arb_midref(acb_realref(z)), 1) > 0 && !arb_contains_si(acb_realref(z), 1);   /* z > 1: t = 1/z on the path */
    if (split) {
        arb_inv(C.c, acb_realref(z), prec);
        arb_sub_ui(C.omc, C.c, 1, prec); arb_neg(C.omc, C.omc);
        C.piece = 1; nb_tanh_sinh(res, nb_nielsen_integrand, &C, prec);
        C.piece = 2; nb_tanh_sinh(part, nb_nielsen_integrand, &C, prec);
        acb_add(res, res, part, prec);
    } else {
        C.piece = 0; nb_tanh_sinh(res, nb_nielsen_integrand, &C, prec);
    }
    acb_rgamma(g, n, prec);                         /* 1/(Gamma(n) Gamma(p + 1)) */
    acb_mul(res, res, g, prec);
    acb_add_ui(part, p, 1, prec);
    acb_rgamma(g, part, prec);
    acb_mul(res, res, g, prec);
    acb_clear(part); acb_clear(g);
    acb_clear(C.n1); acb_clear(C.p); acb_clear(C.z); acb_clear(C.zm1);
    arb_clear(C.c); arb_clear(C.omc);
}

/* The method for one S_{n,p}(z) evaluation (the caller has already decided
 * the call is one the reference evaluates): the Stirling series for a positive
 * integer p and |z| <= 3/4, the Miller series for any other p with |z| <= 1/2
 * (or Re n <= 0, where the integral diverges at t = 0), the quadrature
 * otherwise. */
/* A part whose ball contains zero and whose midpoint is negligible beside the
 * other part is an exact zero: the series' principal powers of positive reals
 * and the quadrature's rounding leave an imaginary dust of ~2^-prec (and a real
 * one on a value that is purely imaginary). */
static void nb_clean_parts(acb_t r, slong prec) {
    arb_ptr part[2] = { acb_realref(r), acb_imagref(r) };
    for (int i = 0; i < 2; i++) {
        arb_ptr a = part[i], b = part[1 - i];
        if (!arb_contains_zero(a) || arb_is_zero(a)) continue;
        mag_t ma, mb; mag_init(ma); mag_init(mb);
        arb_get_mag(ma, a); arb_get_mag(mb, b);
        mag_mul_2exp_si(mb, mb, -(prec - 8));
        if (mag_cmp(ma, mb) <= 0) arb_zero(a);
        mag_clear(ma); mag_clear(mb);
    }
}

static void nb_nielsen_gen(acb_t res, const acb_t n_in, const acb_t p_in, const acb_t z, slong prec) {
    /* local copies: GCC's -Wstringop-overread misreads arb's inline helpers on
     * the const parameters */
    acb_t n, p; acb_init(n); acb_init(p); acb_set(n, n_in); acb_set(p, p_in);
    mag_t mz; mag_init(mz); acb_get_mag(mz, z);
    double az = mag_get_d(mz);
    mag_clear(mz);
    int p_int = acb_is_int(p) && arb_is_positive(acb_realref(p)) &&
                arf_cmpabs_2exp_si(arb_midref(acb_realref(p)), 10) < 0;
    int re_n_pos = arb_is_positive(acb_realref(n));
    if (p_int && az <= 0.75) {
        nb_nielsen(res, n, arf_get_si(arb_midref(acb_realref(p)), ARF_RND_NEAR), z, prec);
    } else if (az < 1.0 && (az <= 0.5 || !re_n_pos)) {
        nb_nielsen_series_gen(res, n, p, z, prec);
    } else if (re_n_pos) {
        nb_nielsen_quad(res, n, p, z, prec);
    } else {
        acb_indeterminate(res);
    }
    if (acb_is_finite(res)) nb_clean_parts(res, prec);
    acb_clear(n); acb_clear(p);
}

Expr* flint_num_nielsen_polylog(const Expr* n, const Expr* p, const Expr* z) {
    return nb_eval3(n, p, z, nb_nielsen_gen);
}

Expr* flint_num_polylog(const Expr* s, const Expr* z) {
    return nb_eval2(s, z, acb_polylog);   /* acb_polylog(w, s, z, prec) */
}

/* Airy Ai/Bi and their derivatives. acb_hypgeom_airy computes all four at once;
 * these adapters pin one output (NULL = skip) so each fits the nb_fn1 shape. The
 * Airy functions are entire -- no branch cuts -- so arb's values match directly. */
static void nb_airy_ai(acb_t r, const acb_t z, slong prec)  { acb_hypgeom_airy(r, NULL, NULL, NULL, z, prec); }
static void nb_airy_aip(acb_t r, const acb_t z, slong prec) { acb_hypgeom_airy(NULL, r, NULL, NULL, z, prec); }
static void nb_airy_bi(acb_t r, const acb_t z, slong prec)  { acb_hypgeom_airy(NULL, NULL, r, NULL, z, prec); }
static void nb_airy_bip(acb_t r, const acb_t z, slong prec) { acb_hypgeom_airy(NULL, NULL, NULL, r, z, prec); }

Expr* flint_num_airy_ai(const Expr* z)       { return nb_eval1(z, nb_airy_ai);  }
/* The error functions are entire, so arb's values need no branch adjustment.
 * Used for large complex |z|, where the MPFR Maclaurin series needs |z|^2/ln2
 * guard bits and ~e|z|^2 terms (16 s at |z| = 100 vs milliseconds here). */
Expr* flint_num_erf(const Expr* z)  { return nb_eval1(z, acb_hypgeom_erf);  }
Expr* flint_num_erfc(const Expr* z) { return nb_eval1(z, acb_hypgeom_erfc); }
Expr* flint_num_barnesg(const Expr* z) { return nb_eval1(z, acb_barnes_g); }
Expr* flint_num_erfi(const Expr* z) { return nb_eval1(z, acb_hypgeom_erfi); }
Expr* flint_num_airy_ai_prime(const Expr* z) { return nb_eval1(z, nb_airy_aip); }
Expr* flint_num_airy_bi(const Expr* z)       { return nb_eval1(z, nb_airy_bi);  }
Expr* flint_num_airy_bi_prime(const Expr* z) { return nb_eval1(z, nb_airy_bip); }

Expr* flint_num_lerchphi(const Expr* z, const Expr* s, const Expr* a) {
    return nb_eval3(z, s, a, acb_dirichlet_lerch_phi);  /* (res, z, s, a, prec) */
}

Expr* flint_num_expint(const Expr* s, const Expr* z) {
    return nb_eval2(s, z, acb_hypgeom_expint);    /* (res, s, z, prec) */
}

static void nb_legp0(acb_t r, const acb_t n, const acb_t m, const acb_t z, slong p) { acb_hypgeom_legendre_p(r, n, m, z, 0, p); }
static void nb_legp1(acb_t r, const acb_t n, const acb_t m, const acb_t z, slong p) { acb_hypgeom_legendre_p(r, n, m, z, 1, p); }
static void nb_legq0(acb_t r, const acb_t n, const acb_t m, const acb_t z, slong p) { acb_hypgeom_legendre_q(r, n, m, z, 0, p); }
static void nb_legq1(acb_t r, const acb_t n, const acb_t m, const acb_t z, slong p) { acb_hypgeom_legendre_q(r, n, m, z, 1, p); }

Expr* flint_num_legendre_p(const Expr* n, const Expr* m, const Expr* z, int type) {
    return nb_eval3(n, m, z, type ? nb_legp1 : nb_legp0);
}
Expr* flint_num_legendre_q(const Expr* n, const Expr* m, const Expr* z, int type) {
    return nb_eval3(n, m, z, type ? nb_legq1 : nb_legq0);
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
/*  Dirichlet L-functions                                              */
/* ------------------------------------------------------------------ */
/* The character arrives as a callback giving chi*(n) = exp(2 Pi I num/den) for
 * a primitive character mod q (src/special_functions/dirichlet.c owns the
 * labelling); it is set on FLINT's OWN generators, so FLINT's Conrey labels
 * never have to agree with Mathematica's. Euler factors (1 - chi*(p) p^-s) for
 * the primes dividing the modulus but not the conductor are applied in arb. */

static int nb_dirichlet_char(dirichlet_group_t G, dirichlet_char_t x, uint64_t q,
                             uint64_t den, flint_chi_exp_fn fn, const void* ctx) {
    dirichlet_group_init(G, (ulong)q);
    dirichlet_char_init(x, G);
    for (slong i = 0; i < G->num; i++) {
        uint64_t num = 0;
        if (!fn(ctx, (uint64_t)G->generators[i], &num)) return 0;
        ulong phi = G->P[i].phi.n;
        /* log = num * phi / den, an integer for a genuine character. */
        mpz_t t; mpz_init_set_ui(t, (unsigned long)num);
        mpz_mul_ui(t, t, (unsigned long)phi);
        int ok = mpz_divisible_ui_p(t, (unsigned long)den);
        if (ok) { mpz_divexact_ui(t, t, (unsigned long)den); x->log[i] = (ulong)mpz_get_ui(t) % phi; }
        mpz_clear(t);
        if (!ok) return 0;
    }
    _dirichlet_char_exp(x, G);
    return 1;
}

/* res *= prod (1 - chi*(p) p^-s). */
static void nb_dirichlet_euler(acb_t res, const acb_t s, int neul, const uint64_t* ep,
                               const uint64_t* enm, uint64_t den, slong prec) {
    acb_t f, c, ps;
    acb_init(f); acb_init(c); acb_init(ps);
    for (int i = 0; i < neul; i++) {
        fmpq_t r; fmpq_init(r);
        fmpq_set_si(r, (slong)(2 * (enm[i] % den)), (ulong)den);
        arb_sin_cos_pi_fmpq(acb_imagref(c), acb_realref(c), r, prec);
        fmpq_clear(r);
        acb_set_ui(ps, (ulong)ep[i]);
        acb_neg(f, s);
        acb_pow(ps, ps, f, prec);
        acb_mul(f, c, ps, prec);
        acb_sub_ui(f, f, 1, prec);
        acb_neg(f, f);
        acb_mul(res, res, f, prec);
    }
    acb_clear(f); acb_clear(c); acb_clear(ps);
}

Expr* flint_num_dirichlet_l(const Expr* s, uint64_t q, uint64_t den,
                            flint_chi_exp_fn fn, const void* ctx,
                            int neul, const uint64_t* eul_p, const uint64_t* eul_num) {
    dirichlet_group_t G; dirichlet_char_t x;
    if (!nb_dirichlet_char(G, x, q, den, fn, ctx)) {
        dirichlet_char_clear(x); dirichlet_group_clear(G);
        return NULL;
    }
    const Expr* args[1] = { s };
    slong outb = pick_out_bits(args, 1);
    acb_t S, R; acb_init(S); acb_init(R);
    Expr* out = NULL;
    for (int attempt = 0; ; attempt++) {
        slong wp = (outb + NB_GUARD) << attempt;
        if (!expr_to_acb(s, S, wp)) break;
        acb_dirichlet_l(R, S, G, x, wp);
        if (neul) nb_dirichlet_euler(R, S, neul, eul_p, eul_num, den, wp);
        if (nb_should_retry(R, outb, attempt)) continue;
        slong bits = nb_emit_bits(R, outb);
        if (bits > 0) out = acb_to_expr(R, bits);
        break;
    }
    acb_clear(S); acb_clear(R);
    dirichlet_char_clear(x); dirichlet_group_clear(G);
    return out;
}

/* Value and first derivative of L(s, chi*) * prod (1 - chi*(p) p^-s) at the
 * ball s: J[0] = value, J[1] = d/ds.
 *
 * For a non-principal chi* (q > 1) the jet is assembled from DEFLATED Hurwitz
 * zeta series, L(s) = q^-s sum_a chi(a) (zeta(s, a/q) - 1/(s-1)), which is
 * exact because sum_a chi(a) = 0 and makes every term entire: a ball containing
 * s = 1 stays finite and wide balls are far less overestimated than with
 * acb_dirichlet_l_jet, whose Hurwitz poles cancel only in exact arithmetic.
 * The principal case (q = 1, Zeta) keeps arb's own jet -- its pole is real.
 * Each Euler factor is the jet (1 - c p^-s, c log(p) p^-s), multiplied in. */
static void nb_dl_jet_hurwitz(acb_struct* J, const acb_t s, const dirichlet_group_t G,
                              const dirichlet_char_t x, slong prec) {
    ulong q = G->q;
    acb_poly_t S, Z; acb_t a, c, acc0, acc1, lq, t;
    acb_poly_init(S); acb_poly_init(Z);
    acb_init(a); acb_init(c); acb_init(acc0); acb_init(acc1); acb_init(lq); acb_init(t);
    acb_poly_set_coeff_acb(S, 0, s);
    acb_poly_set_coeff_si(S, 1, 1);                       /* s + eps */
    for (ulong k = 1; k < q; k++) {
        if (n_gcd(k, q) != 1) continue;
        acb_dirichlet_chi(c, G, x, k, prec);
        acb_set_ui(a, k); acb_div_ui(a, a, q, prec);
        acb_poly_zeta_series(Z, S, a, 1, 2, prec);
        acb_poly_get_coeff_acb(t, Z, 0); acb_addmul(acc0, c, t, prec);
        acb_poly_get_coeff_acb(t, Z, 1); acb_addmul(acc1, c, t, prec);
    }
    /* times q^-s = q^-s (1 - eps log q) */
    acb_set_ui(lq, q); acb_log(lq, lq, prec);
    acb_mul(t, lq, s, prec); acb_neg(t, t); acb_exp(t, t, prec);   /* q^-s */
    acb_mul(J, acc0, t, prec);
    acb_mul(acc1, acc1, t, prec);
    acb_submul(acc1, J, lq, prec);
    acb_set(J + 1, acc1);
    acb_poly_clear(S); acb_poly_clear(Z);
    acb_clear(a); acb_clear(c); acb_clear(acc0); acb_clear(acc1); acb_clear(lq); acb_clear(t);
}

/* The functional equation of a REAL primitive character (root number 1) as a
 * jet: L(s) = (q/Pi)^(1/2-s) Gamma((1-s+d)/2) / Gamma((s+d)/2) L(1-s). Used on
 * the left half-plane, where Hurwitz balls are hopeless. */
static void nb_dl_jet_fe(acb_struct* J, const acb_t s, const dirichlet_group_t G,
                         const dirichlet_char_t x, slong prec) {
    int d = dirichlet_parity_char(G, x);
    acb_poly_t S, A, B, C, P; acb_t t, lqp; acb_struct D[2];
    acb_poly_init(S); acb_poly_init(A); acb_poly_init(B); acb_poly_init(C); acb_poly_init(P);
    acb_init(t); acb_init(lqp); acb_init(D); acb_init(D + 1);
    acb_set_ui(lqp, G->q); acb_const_pi(t, prec); acb_div(lqp, lqp, t, prec); acb_log(lqp, lqp, prec);
    /* A = exp((1/2 - S) log(q/Pi)) */
    acb_one(t); acb_mul_2exp_si(t, t, -1); acb_sub(t, t, s, prec); acb_mul(t, t, lqp, prec);
    acb_poly_set_coeff_acb(A, 0, t);
    acb_neg(t, lqp); acb_poly_set_coeff_acb(A, 1, t);
    acb_poly_exp_series(A, A, 2, prec);
    /* B = Gamma((1 - S + d)/2) */
    acb_set_si(t, 1 + d); acb_sub(t, t, s, prec); acb_mul_2exp_si(t, t, -1);
    acb_poly_set_coeff_acb(B, 0, t);
    acb_set_si(t, -1); acb_mul_2exp_si(t, t, -1); acb_poly_set_coeff_acb(B, 1, t);
    acb_poly_gamma_series(B, B, 2, prec);
    /* C = 1/Gamma((S + d)/2) */
    acb_set_si(t, d); acb_add(t, t, s, prec); acb_mul_2exp_si(t, t, -1);
    acb_poly_set_coeff_acb(C, 0, t);
    acb_one(t); acb_mul_2exp_si(t, t, -1); acb_poly_set_coeff_acb(C, 1, t);
    acb_poly_rgamma_series(C, C, 2, prec);
    /* L(1 - S): value D0, derivative -D1 */
    acb_one(t); acb_sub(t, t, s, prec);
    nb_dl_jet_hurwitz(D, t, G, x, prec);
    acb_poly_fit_length(P, 2);
    acb_poly_set_coeff_acb(P, 0, D);
    acb_neg(t, D + 1); acb_poly_set_coeff_acb(P, 1, t);
    acb_poly_mullow(P, P, A, 2, prec);
    acb_poly_mullow(P, P, B, 2, prec);
    acb_poly_mullow(P, P, C, 2, prec);
    acb_poly_get_coeff_acb(J, P, 0);
    acb_poly_get_coeff_acb(J + 1, P, 1);
    acb_poly_clear(S); acb_poly_clear(A); acb_poly_clear(B); acb_poly_clear(C); acb_poly_clear(P);
    acb_clear(t); acb_clear(lqp); acb_clear(D); acb_clear(D + 1);
}

static void nb_dl_jet(acb_struct* J, const acb_t s, const dirichlet_group_t G,
                      const dirichlet_char_t x, int neul, const uint64_t* ep,
                      const uint64_t* enm, uint64_t den, slong prec) {
    if (G->q == 1) {
        acb_dirichlet_l_jet(J, s, G, x, 0, 2, prec);
    } else {
        arf_t ub; arf_init(ub);
        arb_get_ubound_arf(ub, acb_realref(s), prec);
        /* below -1/4, so that 1 - s stays clear of the s = 1 the Hurwitz form
         * cannot enclose; near 0 the Hurwitz form itself is well behaved. */
        arf_neg(ub, ub);
        int left = arf_cmp_2exp_si(ub, -2) > 0 && dirichlet_order_char(G, x) <= 2;
        arf_clear(ub);
        if (left) nb_dl_jet_fe(J, s, G, x, prec);
        else      nb_dl_jet_hurwitz(J, s, G, x, prec);
    }
    if (!neul) return;
    acb_t c, ps, f0, f1, t, lp;
    acb_init(c); acb_init(ps); acb_init(f0); acb_init(f1); acb_init(t); acb_init(lp);
    for (int i = 0; i < neul; i++) {
        fmpq_t r; fmpq_init(r);
        fmpq_set_si(r, (slong)(2 * (enm[i] % den)), (ulong)den);
        arb_sin_cos_pi_fmpq(acb_imagref(c), acb_realref(c), r, prec);
        fmpq_clear(r);
        acb_set_ui(ps, (ulong)ep[i]);
        acb_log(lp, ps, prec);
        acb_neg(t, s);
        acb_pow(ps, ps, t, prec);                /* p^-s            */
        acb_mul(t, c, ps, prec);                 /* c p^-s          */
        acb_one(f0); acb_sub(f0, f0, t, prec);   /* 1 - c p^-s      */
        acb_mul(f1, t, lp, prec);                /* c log p p^-s    */
        /* (J0, J1) * (f0, f1) */
        acb_mul(t, J + 1, f0, prec);
        acb_addmul(t, J, f1, prec);
        acb_set(J + 1, t);
        acb_mul(J, J, f0, prec);
    }
    acb_clear(c); acb_clear(ps); acb_clear(f0); acb_clear(f1); acb_clear(t); acb_clear(lp);
}

static void nb_dl_widen(const arb_t v, arf_t ylo, arf_t yhi, int* have, slong prec) {
    arf_t t; arf_init(t);
    arb_get_lbound_arf(t, v, prec);
    if (!*have || arf_cmp(t, ylo) < 0) arf_set(ylo, t);
    arb_get_ubound_arf(t, v, prec);
    if (!*have || arf_cmp(t, yhi) > 0) arf_set(yhi, t);
    arf_clear(t);
    *have = 1;
}

/* One certified piece of the enclosure over [lo, hi]; widens [*ylo, *yhi].
 *
 * Arb's L-function jet on a WIDE ball is badly overestimated (it sums Hurwitz
 * zetas whose poles cancel only in exact arithmetic: a ball over [1.1, 3] comes
 * back infinite), so the range is bisected until the derivative ball on a piece
 * excludes 0 -- that piece is monotone and its endpoint values bound it. Near a
 * genuine extremum the pieces never certify; at the depth limit such a piece is
 * bounded by the mean-value form L(mid) +- r sup|L'(piece)|, which is rigorous
 * and, at that width, tight. Bisection only recurses where certification
 * fails, so the cost is O(depth) pieces per extremum. */
static int nb_dl_iv_piece(const arf_t lo, const arf_t hi, int depth,
                          const dirichlet_group_t G, const dirichlet_char_t x,
                          int neul, const uint64_t* ep, const uint64_t* enm, uint64_t den,
                          arf_t ylo, arf_t yhi, int* have, slong prec) {
    acb_t S; acb_struct J[2];
    acb_init(S); acb_init(J); acb_init(J + 1);
    arb_set_interval_arf(acb_realref(S), lo, hi, prec);
    int ok = 0;
    nb_dl_jet(J, S, G, x, neul, ep, enm, den, prec);
    int fin = acb_is_finite(J + 1);
    if (fin && !arb_contains_zero(acb_realref(J + 1))) {
        acb_t P; acb_struct V[2];
        acb_init(P); acb_init(V); acb_init(V + 1);
        ok = 1;
        for (int e = 0; e < 2 && ok; e++) {
            acb_zero(P);
            arb_set_arf(acb_realref(P), e ? hi : lo);
            nb_dl_jet(V, P, G, x, neul, ep, enm, den, prec);
            if (!acb_is_finite(V)) ok = 0;
            else nb_dl_widen(acb_realref(V), ylo, yhi, have, prec);
        }
        acb_clear(P); acb_clear(V); acb_clear(V + 1);
    } else if (depth > 0) {
        arf_t mid; arf_init(mid);
        arf_add(mid, lo, hi, ARF_PREC_EXACT, ARF_RND_DOWN);
        arf_mul_2exp_si(mid, mid, -1);
        ok = nb_dl_iv_piece(lo, mid, depth - 1, G, x, neul, ep, enm, den, ylo, yhi, have, prec)
          && nb_dl_iv_piece(mid, hi, depth - 1, G, x, neul, ep, enm, den, ylo, yhi, have, prec);
        arf_clear(mid);
    } else if (fin) {
        /* L(mid) + [-r, r] * L'(piece) */
        acb_t P; acb_struct V[2]; arb_t rad, t;
        acb_init(P); acb_init(V); acb_init(V + 1); arb_init(rad); arb_init(t);
        acb_zero(P);
        arb_set(acb_realref(P), acb_realref(S));
        mag_zero(arb_radref(acb_realref(P)));               /* the midpoint */
        nb_dl_jet(V, P, G, x, neul, ep, enm, den, prec);
        if (acb_is_finite(V)) {
            arb_zero(rad);
            mag_set(arb_radref(rad), arb_radref(acb_realref(S)));   /* [-r, r] */
            arb_mul(t, rad, acb_realref(J + 1), prec);
            arb_add(t, t, acb_realref(V), prec);
            nb_dl_widen(t, ylo, yhi, have, prec);
            ok = 1;
        }
        acb_clear(P); acb_clear(V); acb_clear(V + 1); arb_clear(rad); arb_clear(t);
    }
    acb_clear(S); acb_clear(J); acb_clear(J + 1);
    return ok;
}

bool flint_dirichlet_l_interval(double lo, double hi, uint64_t q, uint64_t den,
                                flint_chi_exp_fn fn, const void* ctx,
                                int neul, const uint64_t* eul_p, const uint64_t* eul_num,
                                double* out_lo, double* out_hi) {
    if (!(lo <= hi)) return false;
    dirichlet_group_t G; dirichlet_char_t x;
    if (!nb_dirichlet_char(G, x, q, den, fn, ctx)) {
        dirichlet_char_clear(x); dirichlet_group_clear(G);
        return false;
    }
    /* Arb has no finite enclosure of these sums on a ball containing s = 1
     * (the Hurwitz poles cancel only pointwise), so such a range is declined
     * up front -- unevaluated, as Mathematica leaves it -- rather than bisected
     * towards a certain failure. A principal character has a real pole there. */
    if (lo <= 1.0 && 1.0 <= hi) {
        dirichlet_char_clear(x); dirichlet_group_clear(G);
        return false;
    }
    arf_t L, H, ylo, yhi;
    arf_init(L); arf_init(H); arf_init(ylo); arf_init(yhi);
    arf_set_d(L, lo); arf_set_d(H, hi);
    int have = 0;
    int ok = nb_dl_iv_piece(L, H, 12, G, x, neul, eul_p, eul_num, den, ylo, yhi, &have, 96);
    if (ok && have) {
        *out_lo = arf_get_d(ylo, ARF_RND_FLOOR);
        *out_hi = arf_get_d(yhi, ARF_RND_CEIL);
    }
    arf_clear(L); arf_clear(H); arf_clear(ylo); arf_clear(yhi);
    dirichlet_char_clear(x); dirichlet_group_clear(G);
    return ok && have;
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
/*  Tricomi confluent hypergeometric U                                 */
/* ------------------------------------------------------------------ */
/* acb_hypgeom_u has the (res, a, b, z, prec) shape nb_fn3 expects, so this is a
 * straight pass-through: FLINT places the z^(1-b) branch cut on the principal
 * branch (-infinity..0), matching Mathematica, and the precision ladder bounds
 * the whole ball. */
Expr* flint_num_hypergeometric_u(const Expr* a, const Expr* b, const Expr* z) {
    return nb_eval3(a, b, z, acb_hypgeom_u);
}

/* ------------------------------------------------------------------ */
/*  Struve function H_nu(z)                                            */
/* ------------------------------------------------------------------ */
/* FLINT/arb has no Struve kernel, so compose it from the regularized generalized
 * hypergeometric (DLMF 11.2.1):
 *
 *     H_nu(z) = (z/2)^(nu+1) * 1F2~(1; 3/2, nu+3/2; -z^2/4)
 *
 * The regularized 1F2 (acb_hypgeom_pfq with regularized = 1) folds the 1/Gamma
 * factors into the series, so the whole expression is entire in nu -- no pole to
 * dodge at nu = -3/2, -5/2, ... -- and acb_hypgeom_pfq bounds the (convergent)
 * series rigorously. The convergent series loses ~1.44*|z| bits to cancellation;
 * flint_num_struveh folds that allowance into the starting working precision. */
static void nb_struve_h(acb_t res, const acb_t nu, const acb_t z, slong prec) {
    acb_struct a[1], b[2];
    acb_t w, zh, t;
    acb_init(&a[0]); acb_init(&b[0]); acb_init(&b[1]);
    acb_init(w); acb_init(zh); acb_init(t);

    acb_one(&a[0]);                        /* upper parameter {1}            */
    acb_set_d(&b[0], 1.5);                 /* lower parameter 3/2 (exact)    */
    acb_add(&b[1], nu, &b[0], prec);       /* lower parameter nu + 3/2       */

    acb_mul(w, z, z, prec);                /* z^2                            */
    acb_mul_2exp_si(w, w, -2);             /* z^2 / 4                        */
    acb_neg(w, w);                         /* -z^2 / 4                       */

    acb_hypgeom_pfq(res, a, 1, b, 2, w, 1 /* regularized */, prec);

    acb_mul_2exp_si(zh, z, -1);            /* z/2                            */
    acb_add_si(t, nu, 1, prec);            /* nu + 1                         */
    acb_pow(zh, zh, t, prec);              /* (z/2)^(nu+1), principal branch */
    acb_mul(res, res, zh, prec);

    acb_clear(&a[0]); acb_clear(&b[0]); acb_clear(&b[1]);
    acb_clear(w); acb_clear(zh); acb_clear(t);
}

Expr* flint_num_struveh(const Expr* nu, const Expr* z) {
    const Expr* args[2] = { nu, z };
    slong outb = pick_out_bits(args, 2);

    /* Cancellation allowance: read |z| once at a cheap precision and add ~1.45
     * bits per unit of |z| to the starting working precision, so a large |z| is
     * resolved in the first attempt rather than by repeated doubling. This is a
     * precision heuristic only -- the ladder (nb_should_retry / nb_emit_bits)
     * remains the rigorous backstop, so a wrong estimate costs time, not digits. */
    slong cancel = 0;
    {
        acb_t Z0; acb_init(Z0);
        if (expr_to_acb(z, Z0, 53)) {
            arb_t m; arb_init(m);
            acb_abs(m, Z0, 53);
            double za = arf_get_d(arb_midref(m), ARF_RND_UP);
            if (za > 0.0 && za < 1e18) cancel = (slong)(1.45 * za) + 8;
            arb_clear(m);
        }
        acb_clear(Z0);
    }

    acb_t NU, Z, R;
    acb_init(NU); acb_init(Z); acb_init(R);
    Expr* out = NULL;
    for (int attempt = 0; ; attempt++) {
        slong wp = (outb + NB_GUARD + cancel) << attempt;
        if (!expr_to_acb(nu, NU, wp) || !expr_to_acb(z, Z, wp)) break;
        nb_struve_h(R, NU, Z, wp);
        if (nb_should_retry(R, outb, attempt)) continue;
        slong bits = nb_emit_bits(R, outb);
        if (bits > 0) out = acb_to_expr(R, bits);
        break;
    }
    acb_clear(NU); acb_clear(Z); acb_clear(R);
    return out;
}

/* ------------------------------------------------------------------ */
/*  Meijer G-function                                                  */
/* ------------------------------------------------------------------ */
/* FLINT/arb has no Meijer G kernel; it is assembled from Slater's theorem
 * (see src/special_functions/sf_machine_meijerg.h for the formulas):
 *
 *   G = sum_{h<=m} C_h w^{b_h} pFq~(1+b_h-a; 1+b_h-b_{!=h}; (-1)^{p-m-n} w),
 *   C_h = prod_{j<=m,j!=h} pi/sin(pi(b_j-b_h)) prod_{j<=n} Gamma(1+b_h-a_j)
 *                                              / prod_{j>n} Gamma(a_j-b_h),
 *
 * with w = exp(log z / r), the REGULARIZED pFq (acb's `regularized = 1`, so a
 * non-positive-integer lower parameter needs no special case), and the dual
 * G^{n,m}_{q,p}(1-b; 1-a; 1/w) for p > q or p == q, |w| > 1.
 *
 * The log case -- two of b_1..b_m an integer apart, where the pi/sin factors
 * have poles that cancel between terms -- is evaluated by perturbation (the
 * technique mpmath's hypercomb uses): b_j += (j+1) eps with eps = 2^-E,
 * E = target + MG_EPS_GUARD, at a working precision raised by (m-1) E so the
 * cancellation of the O(eps^-(m-1)) poles is absorbed. The perturbation error
 * is O(eps), below the target by construction.
 *
 * The working precision climbs until arb certifies `outb` relative bits (the
 * nb_* ladder rules), with a step sized from the measured shortfall rather than
 * blind doubling, since the dominant loss -- cancellation between Slater terms
 * that grow like exp(c |w|^(1/(q-p))) while G decays -- is predictable from it. */

#define MG_EPS_GUARD 24
#define MG_MAX_ATTEMPTS 16

/* The inputs, either as Exprs (the MPFR / exact-via-N path) or as doubles (the
 * machine fallback); reloaded at each working precision because a Rational is
 * rounded at the precision it is loaded at. */
typedef struct {
    size_t n, p, m, q;
    const Expr* const* a;  const Expr* const* b;  const Expr* z;  const Expr* r;
    const double *ar, *ai, *br, *bi; double r_d, zr, zi;
} MgSrc;

static int mg_load(const MgSrc* s, acb_ptr A, acb_ptr B, acb_t logw, slong wp) {
    acb_t z, r;
    acb_init(z); acb_init(r);
    int ok = 1;
    if (s->a) {
        for (size_t j = 0; j < s->p && ok; j++) ok = expr_to_acb(s->a[j], A + j, wp);
        for (size_t j = 0; j < s->q && ok; j++) ok = expr_to_acb(s->b[j], B + j, wp);
        if (ok) ok = expr_to_acb(s->z, z, wp);
        if (ok) { if (s->r) ok = expr_to_acb(s->r, r, wp); else acb_one(r); }
    } else {
        for (size_t j = 0; j < s->p; j++)
            acb_set_d_d(A + j, s->ar[j], s->ai ? s->ai[j] : 0.0);
        for (size_t j = 0; j < s->q; j++)
            acb_set_d_d(B + j, s->br[j], s->bi ? s->bi[j] : 0.0);
        acb_set_d_d(z, s->zr, s->zi == 0.0 ? 0.0 : s->zi);
        acb_set_d(r, s->r_d);
    }
    if (ok && (acb_is_zero(z) || acb_contains_zero(r))) ok = 0;
    if (ok) {
        acb_log(logw, z, wp);      /* principal: arg(-x) = +pi for exact real z < 0 */
        acb_div(logw, logw, r, wp);
    }
    acb_clear(z); acb_clear(r);
    return ok;
}

/* The regularized pFq of one Slater term, dispatched to arb's specialised
 * kernels where they exist (0F1, 1F1 and 2F1 carry asymptotic expansions and
 * analytic continuation that the generic summation does not). */
static void mg_pfq_reg(acb_t res, acb_srcptr u, slong pu, acb_srcptr l, slong ql,
                       const acb_t x, slong wp) {
    if (pu == 0 && ql == 0) { acb_exp(res, x, wp); return; }
    if (pu == 1 && ql == 0) {               /* (1 - x)^(-u) */
        acb_t t, e; acb_init(t); acb_init(e);
        acb_sub_si(t, x, 1, wp); acb_neg(t, t);
        acb_neg(e, u);
        acb_pow(res, t, e, wp);
        acb_clear(t); acb_clear(e);
        return;
    }
    if (pu == 0 && ql == 1) { acb_hypgeom_0f1(res, l, x, 1, wp); return; }
    if (pu == 1 && ql == 1) { acb_hypgeom_m(res, u, l, x, 1, wp); return; }
    if (pu == 2 && ql == 1) {
        acb_hypgeom_2f1(res, u, u + 1, l, x, ACB_HYPGEOM_2F1_REGULARIZED, wp);
        return;
    }
    acb_hypgeom_pfq(res, u, pu, l, ql, x, 1, wp);
}

/* One Slater sum for G^{m,n}_{p,q}(A; B; w) at log w, working precision wp. */
static void mg_slater(acb_t res, acb_srcptr A, size_t n, size_t p,
                      acb_srcptr B, size_t m, size_t q, const acb_t logw, slong wp) {
    acb_zero(res);
    if (m == 0) return;
    acb_ptr u = _acb_vec_init((slong)(p ? p : 1));
    acb_ptr l = _acb_vec_init((slong)(q ? q : 1));
    acb_t C, t, x, pi, w, sing;
    acb_init(C); acb_init(t); acb_init(x); acb_init(pi); acb_init(w); acb_init(sing);
    acb_const_pi(pi, wp);
    acb_exp(w, logw, wp);
    if (((long)p - (long)m - (long)n) % 2 != 0) acb_neg(x, w); else acb_set(x, w);
    for (size_t h = 0; h < m; h++) {
        acb_one(C);
        for (size_t j = 0; j < m; j++) {
            if (j == h) continue;
            acb_sub(t, B + j, B + h, wp);
            acb_sin_pi(t, t, wp);
            acb_div(C, C, t, wp);
            acb_mul(C, C, pi, wp);
        }
        for (size_t j = 0; j < p; j++) {
            acb_sub(t, B + h, A + j, wp);
            if (j < n) { acb_add_si(t, t, 1, wp); acb_gamma(t, t, wp); }
            else       { acb_neg(t, t); acb_rgamma(t, t, wp); }
            acb_mul(C, C, t, wp);
        }
        for (size_t k = 0; k < p; k++) {
            acb_sub(u + k, B + h, A + k, wp);
            acb_add_si(u + k, u + k, 1, wp);
        }
        size_t nl = 0;
        for (size_t j = 0; j < q; j++) {
            if (j == h) continue;
            acb_sub(l + nl, B + h, B + j, wp);
            acb_add_si(l + nl, l + nl, 1, wp);
            nl++;
        }
        mg_pfq_reg(t, u, (slong)p, l, (slong)nl, x, wp);
        if (!acb_is_finite(t) && p == 2 && nl == 1 && arb_contains_si(acb_realref(x), 1) && arb_contains_zero(acb_imagref(x))) {
            /* x = 1 for a 2F1 term (|w| = 1, p == q == 2): with Re(c-a-b) < 0
             * every term diverges like (1-x)^(c-a-b), but c-a-b = sum a - sum b
             * - 1 is the SAME for every term, so the divergent parts cancel when
             * G is finite there and G(1) is the sum of the regular parts,
             *   2F1~(a,b;c;1-) -> Gamma(c-a-b) / (Gamma(c-a) Gamma(c-b)).
             * The singular coefficients Gamma(a+b-c)/(Gamma(a)Gamma(b)) are
             * summed alongside (into `sing`) and checked to cancel. */
            acb_t s1, s2;
            acb_init(s1); acb_init(s2);
            acb_sub(s1, l, u, wp); acb_sub(s1, s1, u + 1, wp);      /* c-a-b */
            acb_gamma(t, s1, wp);
            acb_sub(s2, l, u, wp); acb_rgamma(s2, s2, wp); acb_mul(t, t, s2, wp);
            acb_sub(s2, l, u + 1, wp); acb_rgamma(s2, s2, wp); acb_mul(t, t, s2, wp);
            acb_neg(s1, s1); acb_gamma(s1, s1, wp);                  /* G(a+b-c) */
            acb_rgamma(s2, u, wp); acb_mul(s1, s1, s2, wp);
            acb_rgamma(s2, u + 1, wp); acb_mul(s1, s1, s2, wp);
            acb_mul(s1, s1, C, wp);
            acb_mul(s2, B + h, logw, wp); acb_exp(s2, s2, wp);
            acb_mul(s1, s1, s2, wp);
            acb_add(sing, sing, s1, wp);
            acb_clear(s1); acb_clear(s2);
        }
        acb_mul(C, C, t, wp);
        acb_mul(t, B + h, logw, wp);
        acb_exp(t, t, wp);
        acb_mul(C, C, t, wp);
        acb_add(res, res, C, wp);
    }
    /* The x = 1 regular-part route is only right when the divergent parts
     * cancel: a sum of singular coefficients that is not (within its own
     * error ball) zero is a genuine singularity of G. */
    if (!acb_is_zero(sing)) {
        arb_t ms, mr; arb_init(ms); arb_init(mr);
        acb_abs(ms, sing, 64); acb_abs(mr, res, 64);
        if (!acb_contains_zero(sing) &&
            !(arf_cmpabs_2exp_si(arb_midref(ms), arf_abs_bound_lt_2exp_si(arb_midref(mr)) - wp / 2) < 0))
            acb_indeterminate(res);
        arb_clear(ms); arb_clear(mr);
    }
    acb_clear(C); acb_clear(t); acb_clear(x); acb_clear(pi); acb_clear(w); acb_clear(sing);
    _acb_vec_clear(u, (slong)(p ? p : 1));
    _acb_vec_clear(l, (slong)(q ? q : 1));
}

/* One evaluation at working precision wp, with the log-case perturbation
 * eps = 2^-E (E ignored when there is no log case). *logc reports whether the
 * log case applies. Returns 0 when the inputs do not load or no expansion
 * applies (p == q >= 3 with |w| = 1). */
static int mg_eval_at(const MgSrc* s, slong wp, slong E, acb_t res, int* logc) {
    size_t p = s->p, q = s->q, mx = (p > q ? p : q) + 1;
    acb_ptr A = _acb_vec_init((slong)mx), B = _acb_vec_init((slong)mx);
    acb_ptr A2 = _acb_vec_init((slong)mx), B2 = _acb_vec_init((slong)mx);
    acb_t logw, t, eps;
    acb_init(logw); acb_init(t); acb_init(eps);
    int ok = 0;
    *logc = 0;
    if (!mg_load(s, A, B, logw, wp)) goto out;
    /* Which expansion: by p vs q, and by |w| when p == q. */
    int dual;
    if (p < q) dual = 0;
    else if (p > q) dual = 1;
    else if (arb_is_negative(acb_realref(logw))) dual = 0;
    else if (arb_is_positive(acb_realref(logw))) dual = 1;
    else {
        if (p > 2) goto out;      /* |w| ~ 1, pF(p-1) with p >= 3: no continuation */
        dual = 0;
    }
    size_t n2, p2, m2, q2;
    if (!dual) {
        n2 = s->n; p2 = p; m2 = s->m; q2 = q;
        for (size_t j = 0; j < p; j++) acb_set(A2 + j, A + j);
        for (size_t j = 0; j < q; j++) acb_set(B2 + j, B + j);
    } else {
        n2 = s->m; p2 = q; m2 = s->n; q2 = p;
        for (size_t j = 0; j < q; j++) { acb_neg(A2 + j, B + j); acb_add_si(A2 + j, A2 + j, 1, wp); }
        for (size_t j = 0; j < p; j++) { acb_neg(B2 + j, A + j); acb_add_si(B2 + j, B2 + j, 1, wp); }
        acb_neg(logw, logw);
    }
    /* Log case: integer-differing b's among the first m2. */
    for (size_t h = 0; h < m2 && !*logc; h++)
        for (size_t j = h + 1; j < m2 && !*logc; j++) {
            acb_sub(t, B2 + j, B2 + h, wp);
            if (acb_contains_int(t)) *logc = 1;
        }
    slong wpe = wp;
    if (*logc) {
        /* (m2 - 1) E absorbs the cancelling eps^-(m2-1) poles of the pi/sin
         * factors; a further 2 E covers a lower parameter 1 + b_h - b_j that was
         * an exact non-positive integer and is now eps away from one (the
         * regularized pFq is continuous there, but arb evaluates it as
         * pFq / Gamma(l) with ~E bits cancelling). */
        wpe = wp + (slong)(m2 + 1) * E;
        acb_one(eps); acb_mul_2exp_si(eps, eps, -E);
        for (size_t j = 0; j < m2; j++) {
            acb_mul_si(t, eps, (slong)(j + 1), wpe);
            acb_add(B2 + j, B2 + j, t, wpe);
        }
    }
    mg_slater(res, A2, n2, p2, B2, m2, q2, logw, wpe);
    ok = 1;
out:
    acb_clear(logw); acb_clear(t); acb_clear(eps);
    _acb_vec_clear(A, (slong)mx); _acb_vec_clear(B, (slong)mx);
    _acb_vec_clear(A2, (slong)mx); _acb_vec_clear(B2, (slong)mx);
    return ok;
}

/* The ladder. Returns 1 with `res` holding a ball certified to `outb` relative
 * bits (or containing zero, or the best reached), 0 when the value does not
 * exist or the inputs are not numeric.
 *
 * In the log case arb's ball does NOT cover the perturbation error, which is
 * O(eps |dG/db|) and can dwarf a small G (MeijerG[{{2,4,31},{3}},{{1.5},
 * {0.1,5}},3] = 1.25e-28 is a sum of terms ~1e-12). So two perturbations are
 * compared, eps and eps 2^-D: their difference estimates the error of the
 * second, and E grows by the measured shortfall until they agree to outb bits. */
static int mg_run(const MgSrc* s, slong outb, acb_t res) {
    acb_t R2, d;
    acb_init(R2); acb_init(d);
    int ok = 0, logc = 0, blowups = 0;
    slong wp = outb + NB_GUARD;
    slong E = outb + MG_EPS_GUARD;
    for (int attempt = 0; attempt < MG_MAX_ATTEMPTS; attempt++) {
        if (!mg_eval_at(s, wp, E, res, &logc)) break;
        if (!acb_is_finite(res)) {
            /* A genuine pole stays one; a ball that merely blew up (Gamma near
             * a pole at low precision) is retried at higher precision. */
            if (++blowups > 2) break;
            wp *= 2;
            continue;
        }
        if (acb_is_zero(res)) { ok = 1; break; }
        /* First the ball: arb's own error bound must reach the target. */
        slong acc = acb_contains_zero(res) ? -wp : acb_rel_accuracy_bits(res);
        if (acc < outb + 4) {
            ok = 1;                      /* best so far, if attempts run out */
            slong step = outb + 4 - acc;
            wp += (step < 32 ? 32 : step) + wp / 4;
            continue;
        }
        ok = 1;
        if (!logc) break;
        /* Then the perturbation (log case): compare eps and eps 2^-D. */
        slong D = outb / 2 + 32;
        if (!mg_eval_at(s, wp, E + D, R2, &logc) || !acb_is_finite(R2)) break;
        slong acc2 = acb_contains_zero(R2) ? -wp : acb_rel_accuracy_bits(R2);
        acb_sub(d, res, R2, wp);
        acb_swap(res, R2);
        if (acc2 < outb + 4) { wp += outb + 4 - acc2 + 32; continue; }
        arb_t ad, ar; arb_init(ad); arb_init(ar);
        acb_abs(ad, d, 64); acb_abs(ar, res, 64);
        double lad = arf_get_d(arb_midref(ad), ARF_RND_UP);
        double lar = arf_get_d(arb_midref(ar), ARF_RND_DOWN);
        arb_clear(ad); arb_clear(ar);
        if (!(lad > 0) || !(lar > 0)) break;          /* identical: converged */
        double agree = log2(lar / lad);               /* bits of the eps result */
        if (agree + (double)D >= (double)outb + 4) break;
        E += (slong)((double)outb + 8 - agree);       /* error is O(eps) */
    }
    acb_clear(R2); acb_clear(d);
    return ok;
}

/* Gauss 2F1(a, b; c; z) at arbitrary precision: acb_hypgeom_2f1, which applies
 * the linear-fractional transformations (z -> 1-z, 1/z, z/(z-1), ...) that a
 * direct series lacks -- the MPFR series HypergeometricPFQ used before needs
 * O(1/(1-|z|)) terms, so 2F1 at |z| = 0.99 cost ~60 ms at 50 digits. Same
 * precision ladder as the other kernels; NULL for non-numeric input, a pole
 * (c a non-positive integer), or a non-finite value. */
Expr* flint_num_hyp2f1(const Expr* a, const Expr* b, const Expr* c, const Expr* z) {
    const Expr* args[4] = { a, b, c, z };
    slong outb = pick_out_bits(args, 4);
    acb_t A, B, C, Z, R;
    acb_init(A); acb_init(B); acb_init(C); acb_init(Z); acb_init(R);
    Expr* out = NULL;
    for (int attempt = 0; ; attempt++) {
        slong wp = (outb + NB_GUARD) << attempt;
        if (!expr_to_acb(a, A, wp) || !expr_to_acb(b, B, wp) ||
            !expr_to_acb(c, C, wp) || !expr_to_acb(z, Z, wp)) break;
        acb_hypgeom_2f1(R, A, B, C, Z, 0, wp);
        if (nb_should_retry(R, outb, attempt)) continue;
        slong bits = nb_emit_bits(R, outb);
        if (bits > 0) out = acb_to_expr(R, bits);
        break;
    }
    acb_clear(A); acb_clear(B); acb_clear(C); acb_clear(Z); acb_clear(R);
    return out;
}

Expr* flint_num_meijerg(const Expr* const* a, size_t n, size_t p,
                        const Expr* const* b, size_t m, size_t q,
                        const Expr* z, const Expr* r) {
    if (n > p || m > q) return NULL;
    long best = numeric_min_inexact_bits(z);
    if (r) { long br = numeric_min_inexact_bits(r); if (br > 0 && (best == 0 || br < best)) best = br; }
    for (size_t j = 0; j < p; j++) {
        long bj = numeric_min_inexact_bits(a[j]);
        if (bj > 0 && (best == 0 || bj < best)) best = bj;
    }
    for (size_t j = 0; j < q; j++) {
        long bj = numeric_min_inexact_bits(b[j]);
        if (bj > 0 && (best == 0 || bj < best)) best = bj;
    }
    if (best < 53) best = 53;
    MgSrc s; memset(&s, 0, sizeof s);
    s.n = n; s.p = p; s.m = m; s.q = q; s.a = a; s.b = b; s.z = z; s.r = r;
    acb_t R; acb_init(R);
    Expr* out = NULL;
    if (mg_run(&s, (slong)best, R)) {
        slong bits = nb_emit_bits(R, (slong)best);
        if (bits > 0) out = acb_to_expr(R, bits);
    }
    acb_clear(R);
    return out;
}

bool flint_mg_eval_d(const double* ar, const double* ai, size_t n, size_t p,
                     const double* br, const double* bi, size_t m, size_t q,
                     double r, double zr, double zi,
                     double* rr, double* ri, bool* is_c) {
    if (n > p || m > q) return false;
    MgSrc s; memset(&s, 0, sizeof s);
    s.n = n; s.p = p; s.m = m; s.q = q;
    s.ar = ar; s.ai = ai; s.br = br; s.bi = bi; s.r_d = r; s.zr = zr; s.zi = zi;
    acb_t R; acb_init(R);
    bool ok = false;
    if (mg_run(&s, 53, R)) {
        *rr = arf_get_d(arb_midref(acb_realref(R)), ARF_RND_NEAR);
        *ri = arf_get_d(arb_midref(acb_imagref(R)), ARF_RND_NEAR);
        *is_c = !arb_is_zero(acb_imagref(R));
        ok = isfinite(*rr) && isfinite(*ri);
    }
    acb_clear(R);
    return ok;
}

/* x_i loaded exactly: (hi + lo) + i (hi + lo) needs at most ~2 * 53 bits plus
 * the exponent gap, well inside the working precision. */
static void mg_laurent_load(acb_t X, const double* x, slong prec) {
    arb_t t; arb_init(t);
    arb_set_d(acb_realref(X), x[0]); arb_set_d(t, x[1]);
    arb_add(acb_realref(X), acb_realref(X), t, prec);
    arb_set_d(acb_imagref(X), x[2]); arb_set_d(t, x[3]);
    arb_add(acb_imagref(X), acb_imagref(X), t, prec);
    arb_clear(t);
}

static int mg_laurent_is_pole(const double* x) {
    return x[1] == 0.0 && x[2] == 0.0 && x[3] == 0.0 && x[0] <= 0.0 && x[0] == floor(x[0]);
}

/* The midpoint of c as a double-double pair. */
static void arb_to_dd(const arb_t c, double* hi, double* lo) {
    arf_t r; arf_init(r);
    *hi = arf_get_d(arb_midref(c), ARF_RND_NEAR);
    arf_set_d(r, *hi);
    arf_sub(r, arb_midref(c), r, ARF_PREC_EXACT, ARF_RND_NEAR);
    *lo = arf_get_d(r, ARF_RND_NEAR);
    arf_clear(r);
}

bool flint_mg_laurent(size_t nf, const double* x, const int* s, const int* e,
                      size_t T, long* v, double* out) {
    if (T == 0) return false;
    /* F = Dn / Nm with Dn the product of 1/Gamma over e = -1 and Nm the
     * product of 1/Gamma over e = +1: both entire in eps, so only their known
     * zero orders (one per pole argument) need dividing out. */
    slong zn = 0, zd = 0;
    for (size_t i = 0; i < nf; i++)
        if (mg_laurent_is_pole(x + 4 * i)) { if (e[i] > 0) zn++; else zd++; }
    slong lenN = (slong)T + zn, lenD = (slong)T + zd;
    bool ok = false;
    for (slong prec = 192; prec <= 768 && !ok; prec *= 2) {
        acb_poly_t h, t, Nm, Dn, F;
        acb_poly_init(h); acb_poly_init(t); acb_poly_init(Nm); acb_poly_init(Dn);
        acb_poly_init(F);
        acb_t X; acb_init(X);
        acb_poly_one(Nm); acb_poly_one(Dn);
        for (size_t i = 0; i < nf; i++) {
            mg_laurent_load(X, x + 4 * i, prec);
            acb_poly_zero(h);
            acb_poly_set_coeff_acb(h, 0, X);
            acb_poly_set_coeff_si(h, 1, s[i] > 0 ? 1 : -1);
            slong len = e[i] > 0 ? lenN : lenD;
            acb_poly_rgamma_series(t, h, len, prec);
            if (mg_laurent_is_pole(x + 4 * i))                /* exact zero at eps = 0 */
                acb_zero(t->coeffs + 0);
            if (e[i] > 0) acb_poly_mullow(Nm, Nm, t, len, prec);
            else          acb_poly_mullow(Dn, Dn, t, len, prec);
        }
        acb_poly_shift_right(Nm, Nm, zn);
        acb_poly_shift_right(Dn, Dn, zd);
        bool good = acb_poly_length(Nm) > 0 && !acb_contains_zero(Nm->coeffs + 0);
        if (good) {
            acb_poly_div_series(F, Dn, Nm, (slong)T, prec);
            /* Every coefficient to >= 110 bits, or negligible beside the largest. */
            arb_t a, big; arb_init(a); arb_init(big);
            arb_zero(big);
            for (slong k = 0; k < (slong)T; k++) {
                acb_t c; acb_init(c);
                acb_poly_get_coeff_acb(c, F, k);
                acb_abs(a, c, 64);
                arb_max(big, big, a, 64);
                acb_clear(c);
            }
            for (slong k = 0; k < (slong)T && good; k++) {
                acb_t c; acb_init(c);
                acb_poly_get_coeff_acb(c, F, k);
                if (!acb_is_finite(c)) good = false;
                else if (!acb_is_exact(c) && acb_rel_accuracy_bits(c) < 110) {
                    mag_t rad, lim; mag_init(rad); mag_init(lim);
                    mag_hypot(rad, arb_radref(acb_realref(c)), arb_radref(acb_imagref(c)));
                    arb_get_mag(lim, big);
                    mag_mul_2exp_si(lim, lim, -115);
                    if (mag_cmp(rad, lim) > 0) good = false;
                    mag_clear(rad); mag_clear(lim);
                }
                if (good) {
                    arb_to_dd(acb_realref(c), &out[4 * k + 0], &out[4 * k + 1]);
                    arb_to_dd(acb_imagref(c), &out[4 * k + 2], &out[4 * k + 3]);
                }
                acb_clear(c);
            }
            arb_clear(a); arb_clear(big);
            ok = good;
        }
        acb_clear(X);
        acb_poly_clear(h); acb_poly_clear(t); acb_poly_clear(Nm); acb_poly_clear(Dn);
        acb_poly_clear(F);
    }
    *v = (long)(zd - zn);
    return ok;
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

/* ------------------------------------------------------------------ */
/*  Certified signs of an exact numeric constant (flint_cert_parts_sign)*/
/* ------------------------------------------------------------------ */

#define CERT_MAX_DEPTH 256

/* Optional variable binding for flint_cert_eval_real (root_trans.c): while set,
 * the symbol g_cert_var is part of the supported language and evaluates to the
 * real ball g_cert_ball, so cert_eval becomes a rigorous INTERVAL evaluator of
 * a univariate expression over that ball.  cert_eval never calls back into the
 * evaluator, so the binding cannot leak into a nested call. */
static const char*      g_cert_var  = NULL;
static const arb_struct* g_cert_ball = NULL;

static const char* cert_head(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return NULL;
    const Expr* h = e->data.function.head;
    return (h && h->type == EXPR_SYMBOL) ? h->data.symbol.name : NULL;
}

static bool cert_is_constant_symbol(const char* s) {
    return s == SYM_Pi || s == SYM_E || s == SYM_EulerGamma || s == SYM_Catalan
        || s == SYM_GoldenRatio || s == SYM_Degree || s == SYM_Khinchin
        || s == SYM_Glaisher;
}

static bool cert_is_int_like(const Expr* e) {
    return e->type == EXPR_INTEGER || e->type == EXPR_BIGINT;
}

static bool cert_is_rational_leaf(const Expr* e) {
    if (cert_is_int_like(e)) return true;
    return cert_head(e) == SYM_Rational && e->data.function.arg_count == 2
        && cert_is_int_like(e->data.function.args[0])
        && cert_is_int_like(e->data.function.args[1]);
}

/* Pure structural pre-pass (no arithmetic): is every node of `e` in the
 * supported constant language?  This is the cheap prefilter — a free symbol or
 * an inexact leaf anywhere is found by pointer compares before any ball is
 * computed, so the hot path (Abs[x], Sign[a b], ...) pays only a short walk. */
static bool cert_supported(const Expr* e, int depth) {
    if (depth > CERT_MAX_DEPTH) return false;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: return true;
        case EXPR_SYMBOL:
            return cert_is_constant_symbol(e->data.symbol.name)
                || (g_cert_var && e->data.symbol.name == g_cert_var);
        case EXPR_FUNCTION: break;
        default: return false;            /* Real / MPFR / String / ... */
    }
    const char* h = cert_head(e);
    if (!h) return false;
    size_t n = e->data.function.arg_count;
    if (h == SYM_Rational) return cert_is_rational_leaf(e);
    if (h == SYM_Complex)
        return n == 2 && cert_is_rational_leaf(e->data.function.args[0])
                      && cert_is_rational_leaf(e->data.function.args[1]);
    if (h == SYM_Plus || h == SYM_Times) { if (n == 0) return false; }
    else if (h == SYM_Power || (h == SYM_Log && n == 2)) { if (n != 2) return false; }
    else if (h == SYM_Sqrt || h == SYM_Exp || h == SYM_Log
          || h == SYM_Sin || h == SYM_Cos || h == SYM_Tan
          || h == SYM_Cot || h == SYM_Sec || h == SYM_Csc
          || h == SYM_Sinh || h == SYM_Cosh || h == SYM_Tanh
          || h == SYM_Coth || h == SYM_Sech || h == SYM_Csch
          || h == SYM_ArcSin || h == SYM_ArcCos || h == SYM_ArcTan
          || h == SYM_ArcSinh || h == SYM_ArcCosh || h == SYM_ArcTanh
          || h == SYM_Gamma || h == SYM_Erf || h == SYM_Erfc
          || h == SYM_Abs || h == SYM_Re || h == SYM_Im || h == SYM_Conjugate) {
        if (n != 1) return false;
    } else {
        return false;
    }
    for (size_t i = 0; i < n; i++)
        if (!cert_supported(e->data.function.args[i], depth + 1)) return false;
    return true;
}

static void cert_set_real(acb_t out, const arb_t r) {
    arb_set(acb_realref(out), r);
    arb_zero(acb_imagref(out));
}

/* Principal Log with a branch-cut guard.  On the cut itself (Im exactly 0,
 * Re certified < 0) the value is Log|z| + I Pi, Mathematica's convention; a
 * ball that merely touches the closed negative real axis is refused. */
static int cert_log(acb_t out, const acb_t z, slong prec) {
    const arb_struct* re = acb_realref(z);
    const arb_struct* im = acb_imagref(z);
    if (arb_is_zero(im)) {
        arb_t t; arb_init(t);
        int ok = 1;
        if (arb_is_positive(re)) {
            arb_log(t, re, prec);
            cert_set_real(out, t);
        } else if (arb_is_negative(re)) {
            arb_neg(t, re);
            arb_log(acb_realref(out), t, prec);
            arb_const_pi(acb_imagref(out), prec);
        } else ok = 0;
        arb_clear(t);
        return ok;
    }
    if (arb_is_positive(re) || !arb_contains_zero(im)) {
        acb_log(out, z, prec);
        return acb_is_finite(out);
    }
    return 0;
}

static int cert_eval(const Expr* e, acb_t out, slong prec);

/* base^x.  Integer x by exact repeated multiplication; a positive real base
 * with a real exponent via arb_pow (real result, exact-zero imaginary part);
 * an exactly-negative real base with exponent p/2 as |b|^(p/2) * I^p (keeps
 * the zero part exactly zero); otherwise Exp[x Log[b]] with the guarded Log. */
static int cert_power(const Expr* e, acb_t out, slong prec) {
    const Expr* be = e->data.function.args[0];
    const Expr* xe = e->data.function.args[1];
    acb_t b, x;
    acb_init(b); acb_init(x);
    int ok = 0;
    if (cert_is_int_like(xe)) {
        if (!cert_eval(be, b, prec)) goto done;
        fmpz_t n; fmpz_init(n);
        expr_to_fmpz(xe, n);
        if (fmpz_sgn(n) < 0 && acb_contains_zero(b)) { fmpz_clear(n); goto done; }
        if (arb_is_zero(acb_imagref(b))) {
            arb_t r; arb_init(r);
            arb_pow_fmpz(r, acb_realref(b), n, prec);
            cert_set_real(out, r);
            arb_clear(r);
        } else {
            acb_pow_fmpz(out, b, n, prec);
        }
        fmpz_clear(n);
        ok = acb_is_finite(out);
        goto done;
    }
    if (be->type == EXPR_SYMBOL && be->data.symbol.name == SYM_E) {
        if (!cert_eval(xe, x, prec)) goto done;
        if (arb_is_zero(acb_imagref(x))) {
            arb_t r; arb_init(r);
            arb_exp(r, acb_realref(x), prec);
            cert_set_real(out, r);
            arb_clear(r);
        } else acb_exp(out, x, prec);
        ok = acb_is_finite(out);
        goto done;
    }
    if (!cert_eval(be, b, prec) || !cert_eval(xe, x, prec)) goto done;
    if (arb_is_zero(acb_imagref(b)) && arb_is_zero(acb_imagref(x))) {
        if (arb_is_positive(acb_realref(b))) {
            arb_t r; arb_init(r);
            arb_pow(r, acb_realref(b), acb_realref(x), prec);
            cert_set_real(out, r);
            arb_clear(r);
            ok = acb_is_finite(out);
            goto done;
        }
        /* Exactly negative base, exponent p/2: |b|^(p/2) * I^p. */
        if (arb_is_negative(acb_realref(b)) && cert_head(xe) == SYM_Rational
            && xe->data.function.args[1]->type == EXPR_INTEGER
            && xe->data.function.args[1]->data.integer == 2
            && xe->data.function.args[0]->type == EXPR_INTEGER) {
            int64_t p = xe->data.function.args[0]->data.integer;
            arb_t r; arb_init(r);
            arb_neg(r, acb_realref(b));
            arb_pow(r, r, acb_realref(x), prec);
            int k = (int)(((p % 4) + 4) % 4);       /* I^p, p odd: I or -I */
            arb_zero(acb_realref(out));
            if (k == 1) arb_set(acb_imagref(out), r);
            else        arb_neg(acb_imagref(out), r);
            arb_clear(r);
            ok = acb_is_finite(out);
            goto done;
        }
    }
    if (acb_contains_zero(b)) goto done;
    {
        acb_t L; acb_init(L);
        if (cert_log(L, b, prec)) {
            acb_mul(L, L, x, prec);
            acb_exp(out, L, prec);
            ok = acb_is_finite(out);
        }
        acb_clear(L);
    }
done:
    acb_clear(b); acb_clear(x);
    return ok;
}

/* One-argument function application f(z).  Real-analytic heads take the arb
 * path on a real ball so the imaginary part stays EXACTLY zero (that exact
 * zero is what certifies realness downstream). */
static int cert_unary(const char* h, acb_t out, const acb_t z, slong prec) {
    bool real = arb_is_zero(acb_imagref(z)) != 0;
    const arb_struct* x = acb_realref(z);
    arb_t r; arb_init(r);
    int ok = 1;
#define CERT_RC(arbf, acbf) do { \
        if (real) { arbf(r, x, prec); cert_set_real(out, r); } \
        else acbf(out, z, prec); } while (0)
    if      (h == SYM_Exp)  CERT_RC(arb_exp,  acb_exp);
    else if (h == SYM_Sin)  CERT_RC(arb_sin,  acb_sin);
    else if (h == SYM_Cos)  CERT_RC(arb_cos,  acb_cos);
    else if (h == SYM_Tan)  CERT_RC(arb_tan,  acb_tan);
    else if (h == SYM_Cot)  CERT_RC(arb_cot,  acb_cot);
    else if (h == SYM_Sec)  CERT_RC(arb_sec,  acb_sec);
    else if (h == SYM_Csc)  CERT_RC(arb_csc,  acb_csc);
    else if (h == SYM_Sinh) CERT_RC(arb_sinh, acb_sinh);
    else if (h == SYM_Cosh) CERT_RC(arb_cosh, acb_cosh);
    else if (h == SYM_Tanh) CERT_RC(arb_tanh, acb_tanh);
    else if (h == SYM_Coth) CERT_RC(arb_coth, acb_coth);
    else if (h == SYM_Sech) CERT_RC(arb_sech, acb_sech);
    else if (h == SYM_Csch) CERT_RC(arb_csch, acb_csch);
    else if (h == SYM_Gamma) CERT_RC(arb_gamma, acb_gamma);
    else if (h == SYM_Erf)  {
        if (real) { arb_hypgeom_erf(r, x, prec); cert_set_real(out, r); }
        else acb_hypgeom_erf(out, z, prec);
    } else if (h == SYM_Erfc) {
        if (real) { arb_hypgeom_erfc(r, x, prec); cert_set_real(out, r); }
        else acb_hypgeom_erfc(out, z, prec);
    } else if (h == SYM_Log || h == SYM_Sqrt) {
        acb_t L; acb_init(L);
        ok = cert_log(L, z, prec);
        if (ok && h == SYM_Sqrt) {
            if (real && arb_is_positive(x)) { arb_sqrt(r, x, prec); cert_set_real(out, r); }
            else if (real) {                       /* exactly negative: I Sqrt[-x] */
                arb_neg(r, x); arb_sqrt(r, r, prec);
                arb_zero(acb_realref(out)); arb_set(acb_imagref(out), r);
            } else acb_sqrt(out, z, prec);        /* off the cut: principal */
        } else if (ok) acb_set(out, L);
        acb_clear(L);
    }
    /* Inverse functions: only where the ball provably avoids the cut. */
    else if (h == SYM_ArcTan || h == SYM_ArcSinh) {    /* cuts on the imaginary axis */
        if (real) {
            if (h == SYM_ArcTan) arb_atan(r, x, prec); else arb_asinh(r, x, prec);
            cert_set_real(out, r);
        } else if (!arb_contains_zero(x)) {
            if (h == SYM_ArcTan) acb_atan(out, z, prec); else acb_asinh(out, z, prec);
        } else ok = 0;
    } else if (h == SYM_ArcSin || h == SYM_ArcCos || h == SYM_ArcTanh || h == SYM_ArcCosh) {
        /* cuts on the real axis: real input strictly inside the domain, or a
         * ball that excludes the real axis altogether. */
        if (real) {
            arb_t a; arb_init(a);
            if (h == SYM_ArcCosh) {
                arb_sub_ui(a, x, 1, prec);
                ok = arb_is_positive(a);
                if (ok) arb_acosh(r, x, prec);
            } else {
                arb_abs(a, x);
                arb_sub_ui(a, a, 1, prec);
                ok = arb_is_negative(a);
                if (ok) {
                    if (h == SYM_ArcSin)      arb_asin(r, x, prec);
                    else if (h == SYM_ArcCos) arb_acos(r, x, prec);
                    else                      arb_atanh(r, x, prec);
                }
            }
            arb_clear(a);
            if (ok) cert_set_real(out, r);
        } else if (!arb_contains_zero(acb_imagref(z))) {
            if (h == SYM_ArcSin)       acb_asin(out, z, prec);
            else if (h == SYM_ArcCos)  acb_acos(out, z, prec);
            else if (h == SYM_ArcTanh) acb_atanh(out, z, prec);
            else                       acb_acosh(out, z, prec);
        } else ok = 0;
    } else if (h == SYM_Abs) {
        acb_abs(r, z, prec); cert_set_real(out, r);
    } else if (h == SYM_Re) {
        cert_set_real(out, acb_realref(z));
    } else if (h == SYM_Im) {
        arb_set(r, acb_imagref(z)); cert_set_real(out, r);
    } else if (h == SYM_Conjugate) {
        acb_conj(out, z);
    } else ok = 0;
#undef CERT_RC
    arb_clear(r);
    return ok && acb_is_finite(out);
}

static int cert_eval(const Expr* e, acb_t out, slong prec) {
    if (e->type == EXPR_SYMBOL) {
        const char* s = e->data.symbol.name;
        if (g_cert_var && s == g_cert_var) {
            cert_set_real(out, g_cert_ball);
            return 1;
        }
        arb_t r; arb_init(r);
        if (s == SYM_Pi) arb_const_pi(r, prec);
        else if (s == SYM_E) arb_const_e(r, prec);
        else if (s == SYM_EulerGamma) arb_const_euler(r, prec);
        else if (s == SYM_Catalan) arb_const_catalan(r, prec);
        else if (s == SYM_Khinchin) arb_const_khinchin(r, prec);
        else if (s == SYM_Glaisher) arb_const_glaisher(r, prec);
        else if (s == SYM_GoldenRatio) {
            arb_sqrt_ui(r, 5, prec); arb_add_ui(r, r, 1, prec); arb_mul_2exp_si(r, r, -1);
        } else if (s == SYM_Degree) {
            arb_const_pi(r, prec); arb_div_ui(r, r, 180, prec);
        } else { arb_clear(r); return 0; }
        cert_set_real(out, r);
        arb_clear(r);
        return 1;
    }
    if (cert_is_rational_leaf(e)) {
        arb_zero(acb_imagref(out));
        return scalar_to_arb(e, acb_realref(out), prec);
    }
    const char* h = cert_head(e);
    if (!h) return 0;
    size_t n = e->data.function.arg_count;
    if (h == SYM_Complex)
        return scalar_to_arb(e->data.function.args[0], acb_realref(out), prec)
            && scalar_to_arb(e->data.function.args[1], acb_imagref(out), prec);
    if (h == SYM_Plus || h == SYM_Times) {
        acb_t t; acb_init(t);
        int ok = cert_eval(e->data.function.args[0], out, prec);
        for (size_t i = 1; ok && i < n; i++) {
            ok = cert_eval(e->data.function.args[i], t, prec);
            if (ok) {
                if (h == SYM_Plus) acb_add(out, out, t, prec);
                else               acb_mul(out, out, t, prec);
            }
        }
        acb_clear(t);
        return ok && acb_is_finite(out);
    }
    if (h == SYM_Power) return cert_power(e, out, prec);
    if (h == SYM_Log && n == 2) {                  /* Log[b, x] = Log[x]/Log[b] */
        acb_t b, x, lb;
        acb_init(b); acb_init(x); acb_init(lb);
        int ok = cert_eval(e->data.function.args[0], b, prec)
              && cert_eval(e->data.function.args[1], x, prec)
              && cert_log(lb, b, prec) && cert_log(out, x, prec)
              && !acb_contains_zero(lb);
        if (ok) {
            if (arb_is_zero(acb_imagref(lb)) && arb_is_zero(acb_imagref(out))) {
                arb_div(acb_realref(out), acb_realref(out), acb_realref(lb), prec);
            } else acb_div(out, out, lb, prec);
            ok = acb_is_finite(out);
        }
        acb_clear(b); acb_clear(x); acb_clear(lb);
        return ok;
    }
    acb_t z; acb_init(z);
    int ok = cert_eval(e->data.function.args[0], z, prec)
          && cert_unary(h, out, z, prec);
    acb_clear(z);
    return ok;
}

static int cert_arb_sign(const arb_t x) {
    if (arb_is_zero(x)) return 0;
    if (arb_is_positive(x)) return 1;
    if (arb_is_negative(x)) return -1;
    return 2;
}

bool flint_cert_parts_sign(const Expr* e, int* re_sgn, int* im_sgn) {
    *re_sgn = *im_sgn = 2;
    if (!e || !cert_supported(e, 0)) return false;
    static const slong ladder[] = { 64, 256, 1024 };
    bool any = false;
    acb_t z; acb_init(z);
    for (size_t i = 0; i < sizeof ladder / sizeof ladder[0]; i++) {
        if (!cert_eval(e, z, ladder[i])) continue;   /* a wider ball may refuse a cut */
        any = true;
        *re_sgn = cert_arb_sign(acb_realref(z));
        *im_sgn = cert_arb_sign(acb_imagref(z));
        if (*re_sgn != 2 && *im_sgn != 2) break;
    }
    acb_clear(z);
    return any;
}

void flint_cert_binding_clear(void) { g_cert_var = NULL; g_cert_ball = NULL; }

int flint_cert_eval_real(const Expr* e, const char* var, const arb_t xball,
                         arb_t out, slong prec) {
    if (!e) return 0;
    const char* sv = g_cert_var; const arb_struct* sb = g_cert_ball;
    g_cert_var = var; g_cert_ball = xball;
    int ok = 0;
    if (cert_supported(e, 0)) {
        acb_t z; acb_init(z);
        if (cert_eval(e, z, prec) && arb_is_zero(acb_imagref(z))
            && arb_is_finite(acb_realref(z))) {
            arb_set(out, acb_realref(z));
            ok = 1;
        }
        acb_clear(z);
    }
    g_cert_var = sv; g_cert_ball = sb;
    return ok;
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

void flint_cert_binding_clear(void) {}

bool flint_cert_parts_sign(const Expr* e, int* re_sgn, int* im_sgn) {
    (void)e; *re_sgn = *im_sgn = 2; return false;
}

Expr* flint_num_zeta(const Expr* s) { (void)s; return NULL; }
Expr* flint_num_hurwitz_zeta(const Expr* s, const Expr* a) { (void)s; (void)a; return NULL; }
Expr* flint_num_polygamma(const Expr* n, const Expr* z) { (void)n; (void)z; return NULL; }
Expr* flint_num_polygamma_frac(const Expr* nu, const Expr* z) { (void)nu; (void)z; return NULL; }
Expr* flint_num_polylog(const Expr* s, const Expr* z) { (void)s; (void)z; return NULL; }
Expr* flint_num_airy_ai(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_erf(const Expr* z)  { (void)z; return NULL; }
Expr* flint_num_erfc(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_barnesg(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_erfi(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_airy_ai_prime(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_airy_bi(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_airy_bi_prime(const Expr* z) { (void)z; return NULL; }
Expr* flint_num_lerchphi(const Expr* z, const Expr* s, const Expr* a) { (void)z; (void)s; (void)a; return NULL; }
Expr* flint_num_stieltjes(const Expr* n, const Expr* a) { (void)n; (void)a; return NULL; }
Expr* flint_num_expint(const Expr* s, const Expr* z) { (void)s; (void)z; return NULL; }
Expr* flint_num_nielsen_polylog(const Expr* n, const Expr* p, const Expr* z) { (void)n; (void)p; (void)z; return NULL; }
Expr* flint_num_bessel_j_zero(const Expr* nu, double seed, long bits) { (void)nu; (void)seed; (void)bits; return NULL; }
Expr* flint_num_qpochhammer(const Expr* a, const Expr* q, const Expr* n) { (void)a; (void)q; (void)n; return NULL; }
Expr* flint_num_legendre_p(const Expr* n, const Expr* m, const Expr* z, int type) { (void)n; (void)m; (void)z; (void)type; return NULL; }
Expr* flint_num_legendre_q(const Expr* n, const Expr* m, const Expr* z, int type) { (void)n; (void)m; (void)z; (void)type; return NULL; }
Expr* flint_num_elliptic_k(const Expr* m) { (void)m; return NULL; }
Expr* flint_num_elliptic_e(const Expr* m) { (void)m; return NULL; }
Expr* flint_num_elliptic_f(const Expr* phi, const Expr* m) { (void)phi; (void)m; return NULL; }
Expr* flint_num_elliptic_e_inc(const Expr* phi, const Expr* m) { (void)phi; (void)m; return NULL; }
Expr* flint_num_elliptic_pi(const Expr* n, const Expr* m) { (void)n; (void)m; return NULL; }
Expr* flint_num_elliptic_pi_inc(const Expr* n, const Expr* phi, const Expr* m) {
    (void)n; (void)phi; (void)m; return NULL;
}
Expr* flint_num_hypergeometric_u(const Expr* a, const Expr* b, const Expr* z) {
    (void)a; (void)b; (void)z; return NULL;
}
Expr* flint_num_struveh(const Expr* nu, const Expr* z) { (void)nu; (void)z; return NULL; }
Expr* flint_num_hyp2f1(const Expr* a, const Expr* b, const Expr* c, const Expr* z) {
    (void)a; (void)b; (void)c; (void)z; return NULL;
}
Expr* flint_num_meijerg(const Expr* const* a, size_t n, size_t p,
                        const Expr* const* b, size_t m, size_t q,
                        const Expr* z, const Expr* r) {
    (void)a; (void)n; (void)p; (void)b; (void)m; (void)q; (void)z; (void)r; return NULL;
}
bool flint_mg_eval_d(const double* ar, const double* ai, size_t n, size_t p,
                     const double* br, const double* bi, size_t m, size_t q,
                     double r, double zr, double zi,
                     double* rr, double* ri, bool* is_c) {
    (void)ar; (void)ai; (void)n; (void)p; (void)br; (void)bi; (void)m; (void)q;
    (void)r; (void)zr; (void)zi; (void)rr; (void)ri; (void)is_c; return false;
}
bool flint_mg_laurent(size_t nf, const double* x, const int* s, const int* e,
                      size_t T, long* v, double* out) {
    (void)nf; (void)x; (void)s; (void)e; (void)T; (void)v; (void)out; return false;
}
Expr* flint_num_dirichlet_l(const Expr* s, uint64_t q, uint64_t den,
                            flint_chi_exp_fn fn, const void* ctx,
                            int neul, const uint64_t* eul_p, const uint64_t* eul_num) {
    (void)s; (void)q; (void)den; (void)fn; (void)ctx; (void)neul; (void)eul_p; (void)eul_num;
    return NULL;
}
bool flint_dirichlet_l_interval(double lo, double hi, uint64_t q, uint64_t den,
                                flint_chi_exp_fn fn, const void* ctx,
                                int neul, const uint64_t* eul_p, const uint64_t* eul_num,
                                double* out_lo, double* out_hi) {
    (void)lo; (void)hi; (void)q; (void)den; (void)fn; (void)ctx; (void)neul; (void)eul_p;
    (void)eul_num; (void)out_lo; (void)out_hi; return false;
}
void  flint_num_bridge_init(void) { /* no FLINT/MPFR: nothing to register */ }

#endif

/* ------------------------------------------------------------------ */
/*  The numeric tail of a builtin (always compiled)                    */
/* ------------------------------------------------------------------ */

static bool nb_is_number(const Expr* e) {
    if (e->type == EXPR_INTEGER || e->type == EXPR_BIGINT || e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL ||
        e->data.function.arg_count != 2) return false;
    const char* h = e->data.function.head->data.symbol.name;
    if (h == SYM_Rational) return true;
    return h == SYM_Complex && nb_is_number(e->data.function.args[0]) &&
           nb_is_number(e->data.function.args[1]);
}

bool flint_num_call_inexact(const Expr* call) {
    if (!call || call->type != EXPR_FUNCTION) return false;
    bool inexact = false;
    for (size_t i = 0; i < call->data.function.arg_count; i++) {
        const Expr* a = call->data.function.args[i];
        if (!nb_is_number(a)) return false;
        if (arg_is_inexact(a)) inexact = true;
    }
    return inexact;
}

Expr* flint_num_for_call(Expr* r, const Expr* call) {
    if (!r) return NULL;
#ifdef USE_MPFR
    long bits = 0;
    for (size_t i = 0; i < call->data.function.arg_count; i++) {
        long b = numeric_min_inexact_bits(call->data.function.args[i]);
        if (b > 0 && (bits == 0 || b < bits)) bits = b;
    }
    return bits <= 53 ? numeric_lower_mpfr_to_machine(r) : r;
#else
    (void)call;
    return r;
#endif
}
