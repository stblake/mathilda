/* Mathilda -- WhittakerM, the Whittaker M function M_{k,m}(z).
 *
 *   WhittakerM[k, m, z] = e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z].
 *
 * Evaluation is layered so each argument shape takes the cheapest numeric route
 * (mirroring the LaguerreL / HermiteH confluent-hypergeometric design):
 *
 *   z == 0 (exact), Re(m) > -1/2   ->  exact 0 (M_{k,m} ~ z^(m+1/2) near 0).
 *   machine-precision call         ->  the shared double helpers
 *   (k, m exact or machine real;         sf_machine_whittaker_m (real z > 0) /
 *    z machine real or complex)          sf_machine_whittaker_m_c (complex z, or
 *                                        real z < 0 on the principal branch),
 *                                      bit-for-bit with the NDArray element-wise
 *                                      buffer loop and the Compile[] VM.
 *                                      Declines (2m+1 pole, cancellation,
 *                                       overflow) -> template.
 *   some inexact argument          ->  the confluent-hypergeometric form built as
 *                                        e^(-z/2) z^(m+1/2) Hypergeometric1F1[...],
 *                                      evaluated so it inherits 1F1's machine +
 *                                      MPFR, real + complex numerics. Kept only
 *                                      if it collapses to a number.
 *   everything else (symbolic)     ->  stays symbolic (return NULL).
 *
 * WhittakerM is a single growing PRODUCT (prefactor * 1F1), not a difference, so
 * it adds no cancellation over 1F1 itself and needs no guard digits (unlike
 * HermiteH): bad 1F1 cancellation is declined inside sf_machine_pfq / re-summed
 * in MPFR by the 1F1 evaluator.
 *
 * Attributes: Listable, NumericFunction, Protected.
 *
 * Deferred / best-effort, as the sibling special functions: tight Interval
 * enclosures (the generic NumericFunction + D[] interval threading in eval.c
 * declines here because D[WhittakerM,z] mentions WhittakerM itself -- the
 * self-referential-derivative guard), the general asymptotic Series at infinity
 * (a two-scale Stokes transseries; only the terminating/elementary case is
 * emitted, in calculus/series.c), symbolic FunctionExpand, and general-n
 * SeriesCoefficient. The D[] rule and the origin Series live in
 * src/calculus/deriv.c and src/calculus/series.c respectively.
 */
#include "whittaker.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <gmp.h>

#include "arithmetic.h"      /* is_complex */
#include "attr.h"
#include "common.h"          /* builtin_arg_error */
#include "eval.h"            /* eval_and_free */
#include "expr.h"
#include "ndarray.h"         /* is_ndarray, ndt_get/_set, expr_new_ndarray_like, ... */
#include "sf_machine.h"      /* sf_machine_whittaker_m */
#include "sf_machine_hyper.h" /* complex twins */
#include "sf_scalar.h"       /* sf_machine_scalar_args, sf_machine_leaf */
#include "symtab.h"
#include "../sym_names.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/* Expression-building shorthands                                      */
/* ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}
/* 1/2 and -1/2 as exact rationals. */
static Expr* mk_half(void)     { return mk_fn2("Rational", expr_new_integer(1),  expr_new_integer(2)); }
static Expr* mk_neg_half(void) { return mk_fn2("Rational", expr_new_integer(-1), expr_new_integer(2)); }

/* ------------------------------------------------------------------ */
/* Argument-shape predicates / extractors                             */
/* ------------------------------------------------------------------ */

/* True if `e` is an inexact numeric leaf (Real or MPFR). */
static bool whit_is_inexact_leaf(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
#ifdef USE_MPFR
    if (e->type == EXPR_MPFR) return true;
#endif
    return false;
}

/* True if `e` carries an inexact numeric leaf (also inside Complex[..]). */
static bool whit_is_inexact(const Expr* e) {
    if (whit_is_inexact_leaf(e)) return true;
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im))
        return whit_is_inexact_leaf(re) || whit_is_inexact_leaf(im);
    return false;
}

/* True iff `e` is a pure numeric quantity (Integer/Real/MPFR/BigInt, or
 * Rational/Complex of such) -- i.e. the result actually collapsed to a number. */
static bool whit_is_number(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER:
        case EXPR_REAL:
        case EXPR_BIGINT:
#ifdef USE_MPFR
        case EXPR_MPFR:
#endif
            return true;
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        (e->data.function.head->data.symbol.name == SYM_Rational ||
         e->data.function.head->data.symbol.name == SYM_Complex)) {
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            if (!whit_is_number(e->data.function.args[i])) return false;
        return true;
    }
    return false;
}

/* Value of a real numeric leaf as a double (Integer/Real/BigInt/MPFR/Rational). */
static bool whit_real_number_double(const Expr* e, double* out) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: *out = (double)e->data.integer; return true;
        case EXPR_REAL:    *out = e->data.real;            return true;
        case EXPR_BIGINT:  *out = mpz_get_d(e->data.bigint); return true;
#ifdef USE_MPFR
        case EXPR_MPFR:    *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true;
#endif
        default: break;
    }
    if (e->type == EXPR_FUNCTION && e->data.function.head &&
        e->data.function.head->type == EXPR_SYMBOL &&
        e->data.function.head->data.symbol.name == SYM_Rational &&
        e->data.function.arg_count == 2) {
        double p, q;
        if (whit_real_number_double(e->data.function.args[0], &p) &&
            whit_real_number_double(e->data.function.args[1], &q) && q != 0.0) {
            *out = p / q;
            return true;
        }
    }
    return false;
}

/* Real part of a numeric `e` as a double (for the z == 0 special value).
 * Declines a symbolic or non-numeric argument. */
static bool whit_real_part_double(const Expr* e, double* out) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) return whit_real_number_double(re, out);
    return whit_real_number_double(e, out);
}

/* ------------------------------------------------------------------ */
/* Confluent-hypergeometric template                                   */
/* ------------------------------------------------------------------ */

/* Build e^(-z/2) z^(m+1/2) Hypergeometric1F1[m-k+1/2, 2m+1, z] as an expression
 * and evaluate it, so it inherits 1F1's machine + MPFR, real + complex numerics
 * (and returns a Complex result for z < 0 or complex parameters). Returns NULL
 * if the result did not collapse to a number (symbolic k / m). */
static Expr* whittaker_m_numeric_hyp(const Expr* k_e, const Expr* m_e, const Expr* z_e) {
    /* a = m - k + 1/2 */
    Expr* a = mk_fn3("Plus",
                     expr_copy((Expr*)m_e),
                     mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)k_e)),
                     mk_half());
    /* b = 2 m + 1 */
    Expr* b = mk_fn2("Plus",
                     mk_fn2("Times", expr_new_integer(2), expr_copy((Expr*)m_e)),
                     expr_new_integer(1));
    Expr* hyp = mk_fn3("Hypergeometric1F1", a, b, expr_copy((Expr*)z_e));

    /* e^(-z/2) */
    Expr* expfac = mk_fn1("Exp", mk_fn2("Times", mk_neg_half(), expr_copy((Expr*)z_e)));
    /* z^(m+1/2) */
    Expr* zpow = mk_fn2("Power", expr_copy((Expr*)z_e),
                        mk_fn2("Plus", expr_copy((Expr*)m_e), mk_half()));

    Expr* prod = mk_fn3("Times", expfac, zpow, hyp);
    Expr* out = eval_and_free(prod);
    if (out && whit_is_number(out)) return out;
    expr_free(out);
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Machine-precision fast path (shared by the scalar and buffer paths) */
/* ------------------------------------------------------------------ */

typedef bool (*WhitReal)(double, double, double, double*);
typedef bool (*WhitCplx)(double, double, double, double, double*, double*);

/* One point: real k, m; z real or complex. A real z <= 0 is on the branch cut
 * of z^(m+1/2): the real kernel declines and the complex one answers on the
 * principal branch (arg z = +pi), as Mathematica does. */
static bool whit_point(WhitReal fr, WhitCplx fc, double k, double m,
                       double zr, double zi, bool zc,
                       double* rr, double* ri, bool* is_c) {
    if (!zc) {
        *is_c = false; *ri = 0.0;
        if (fr(k, m, zr, rr)) return true;
        if (!(zr < 0.0)) return false;
        zi = 0.0;
    }
    *is_c = true;
    return fc(k, m, zr, zi, rr, ri);
}

/* The scalar gate: k, m, z machine-scalar leaves (Integer / Rational / machine
 * Real, z possibly Complex), at least one machine Real, k and m real. Shares
 * the kernels with the buffer path below and the Compile[] VM. */
static Expr* whit_machine(Expr* res, WhitReal fr, WhitCplx fc) {
    SFMachineArgs ma;
    if (!sf_machine_scalar_args(res, &ma) || ma.n != 3) return NULL;
    if (ma.is_complex[0] || ma.is_complex[1]) return NULL;
    double rr, ri; bool is_c;
    if (!whit_point(fr, fc, ma.re[0], ma.re[1], ma.re[2], ma.im[2], ma.is_complex[2],
                    &rr, &ri, &is_c)) return NULL;
    return is_c ? sf_machine_complex_result(rr, ri) : expr_new_real(rr);
}

/* The buffer path: thread the same kernels over a visible NDArray / packed List
 * z, k and m held constant (exact or machine). Once z IS an NDArray this must
 * not return NULL -- the scalar driver cannot consume a buffer -- so anything it
 * cannot do (an int64 buffer, an MPFR or complex k/m, a declined element) is
 * handed to ndarray_delist_and_reeval, where each element takes the accurate
 * scalar route. */
static Expr* whit_ndarray(WhitReal fr, WhitCplx fc, const Expr* k_e, const Expr* m_e,
                          const Expr* z_e, Expr* res) {
    if (!is_ndarray(z_e)) return NULL;             /* scalar driver handles */
    NDType dt = z_e->data.ndarray.dtype;
    double kd, md; bool ixk = false, ixm = false;
    if (dt == NDT_INT64 || dt == NDT_BOOL ||
        !sf_machine_leaf(k_e, &kd, &ixk) || !sf_machine_leaf(m_e, &md, &ixm))
        return ndarray_delist_and_reeval(res);
    bool zc = ndt_is_complex(dt);
    size_t sz = ndarray_size(z_e);
    if (dt == NDT_FLOAT64) {                       /* the common case: real in, real out */
        const double* x = (const double*)z_e->data.ndarray.data;
        double* o = malloc(sizeof(double) * (sz ? sz : 1));
        if (!o) return ndarray_delist_and_reeval(res);
        size_t i = 0;
        for (; i < sz; i++) if (!fr(kd, md, x[i], &o[i])) break;
        if (i == sz)
            return expr_new_ndarray_like(z_e, z_e->data.ndarray.rank,
                                         z_e->data.ndarray.dims, o, NDT_FLOAT64);
        free(o);                                   /* a cut / declined point: general loop */
    }
    double* buf = malloc(sizeof(double) * 2 * (sz ? sz : 1));
    if (!buf) return ndarray_delist_and_reeval(res);
    bool any_c = zc;
    for (size_t i = 0; i < sz; i++) {
        double xr, xi, rr, ri; bool is_c;
        ndt_get(z_e->data.ndarray.data, i, dt, &xr, &xi);
        if (!whit_point(fr, fc, kd, md, xr, xi, zc, &rr, &ri, &is_c)) {
            free(buf);
            return ndarray_delist_and_reeval(res);
        }
        buf[2 * i] = rr; buf[2 * i + 1] = ri;
        any_c |= is_c;
    }
    NDType odt = any_c ? (dt == NDT_FLOAT32 || dt == NDT_COMPLEX32 ? NDT_COMPLEX32
                                                                    : NDT_COMPLEX64)
                       : dt;
    void* out = malloc(ndt_elem_size(odt) * (sz ? sz : 1));
    if (!out) { free(buf); return ndarray_delist_and_reeval(res); }
    for (size_t i = 0; i < sz; i++) ndt_set(out, i, odt, buf[2 * i], buf[2 * i + 1]);
    free(buf);
    return expr_new_ndarray_like(z_e, z_e->data.ndarray.rank, z_e->data.ndarray.dims,
                                 out, odt);
}

/* ------------------------------------------------------------------ */
/* Scalar driver                                                       */
/* ------------------------------------------------------------------ */

static Expr* whittaker_m_eval(Expr* k_e, Expr* m_e, Expr* z_e) {
    /* Special value: M_{k,m}(0) = 0 for Re(m) > -1/2 (leading power z^(m+1/2)).
     * Only the exact integer 0; a Real 0.0 flows through the template to 0.0. */
    if (z_e && z_e->type == EXPR_INTEGER && z_e->data.integer == 0) {
        double rem;
        if (whit_real_part_double(m_e, &rem) && rem > -0.5)
            return expr_new_integer(0);
        /* symbolic m, or Re(m) <= -1/2 (singular): fall through / stay symbolic */
    }

    /* Numericise only when some argument is inexact; otherwise stay symbolic
     * (so WhittakerM[2,3,x] and all-exact WhittakerM[2,3,4] do not expand). */
    if (!(whit_is_inexact(k_e) || whit_is_inexact(m_e) || whit_is_inexact(z_e)))
        return NULL;

    /* General confluent-hypergeometric template (MPFR precision, complex args,
     * z <= 0, exact-rational k/m with an inexact z). */
    Expr* out = whittaker_m_numeric_hyp(k_e, m_e, z_e);
    if (out) return out;
    return NULL;                               /* stay symbolic */
}

/* ================================================================== */
/* WhittakerW: W_{k,m}(z) = e^(-z/2) z^(m+1/2) HypergeometricU[m-k+1/2, 2m+1, z]. */
/* ================================================================== */

/* True iff `e` is a real number (not a Complex[..] with a non-zero imaginary
 * part). A bare Integer/Real/Rational is real; Complex[re, 0] is real too. */
static bool whit_is_real_number(const Expr* e) {
    Expr *re, *im;
    if (is_complex((Expr*)e, &re, &im)) {
        double imd;
        return whit_real_number_double(im, &imd) && imd == 0.0;
    }
    return true;
}

/* Build e^(-z/2) z^(m+1/2) HypergeometricU[m-k+1/2, 2m+1, z] and evaluate it, so
 * WhittakerW inherits HypergeometricU's machine + MPFR, real + complex numerics
 * (exactly as whittaker_m_numeric_hyp rides Hypergeometric1F1). Returns NULL if
 * the result did not collapse to a number. */
static Expr* whittaker_w_numeric_hyp(const Expr* k_e, const Expr* m_e, const Expr* z_e) {
    Expr* a = mk_fn3("Plus",
                     expr_copy((Expr*)m_e),
                     mk_fn2("Times", expr_new_integer(-1), expr_copy((Expr*)k_e)),
                     mk_half());                                  /* m - k + 1/2 */
    Expr* b = mk_fn2("Plus",
                     mk_fn2("Times", expr_new_integer(2), expr_copy((Expr*)m_e)),
                     expr_new_integer(1));                        /* 2 m + 1 */
    Expr* u = mk_fn3("HypergeometricU", a, b, expr_copy((Expr*)z_e));
    Expr* expfac = mk_fn1("Exp", mk_fn2("Times", mk_neg_half(), expr_copy((Expr*)z_e)));
    Expr* zpow = mk_fn2("Power", expr_copy((Expr*)z_e),
                        mk_fn2("Plus", expr_copy((Expr*)m_e), mk_half()));
    Expr* prod = mk_fn3("Times", expfac, zpow, u);
    Expr* out = eval_and_free(prod);
    if (out && whit_is_number(out)) return out;
    expr_free(out);
    return NULL;
}

static Expr* whittaker_w_eval(Expr* k_e, Expr* m_e, Expr* z_e) {
    /* Special values at z == 0 (exact Integer 0 only; a Real 0.0 flows through).
     * W near 0 has two Puiseux branches z^(1/2 +- m), so the leading exponent is
     * 1/2 - |Re m|: for |Re m| < 1/2 it vanishes, for |Re m| > 1/2 it blows up,
     * and at |Re m| = 1/2 (m = +-1/2, real) it is the finite 1/Gamma[1-k]. */
    if (z_e && z_e->type == EXPR_INTEGER && z_e->data.integer == 0) {
        double rem;
        if (whit_real_part_double(m_e, &rem)) {
            double am = fabs(rem);
            if (am < 0.5) return expr_new_integer(0);
            if (am > 0.5) return expr_new_symbol("ComplexInfinity");
            /* |Re m| == 1/2: finite 1/Gamma[1-k] only for a real m = +-1/2. */
            if (whit_is_real_number(m_e))
                return eval_and_free(
                    mk_fn2("Power",
                           mk_fn1("Gamma",
                                  mk_fn2("Plus", expr_new_integer(1),
                                         mk_fn2("Times", expr_new_integer(-1),
                                                expr_copy(k_e)))),
                           expr_new_integer(-1)));
            /* complex m on the boundary: no limit -> stay symbolic */
        }
        /* symbolic m: stay symbolic */
    }

    /* Numericise only when some argument is inexact; otherwise stay symbolic. */
    if (!(whit_is_inexact(k_e) || whit_is_inexact(m_e) || whit_is_inexact(z_e)))
        return NULL;

    Expr* out = whittaker_w_numeric_hyp(k_e, m_e, z_e);
    if (out) return out;
    return NULL;                               /* stay symbolic */
}

Expr* builtin_whittaker_w(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc != 3)
        return builtin_arg_error("WhittakerW", argc, 3, 3);
    Expr* nd = whit_ndarray(sf_machine_whittaker_w, sf_machine_whittaker_w_c,
                            args[0], args[1], args[2], res);
    if (nd) return nd;
    Expr* m = whit_machine(res, sf_machine_whittaker_w, sf_machine_whittaker_w_c);
    if (m) return m;
    return whittaker_w_eval(args[0], args[1], args[2]);
}

/* ------------------------------------------------------------------ */
/* Builtin entry point                                                 */
/* ------------------------------------------------------------------ */

Expr* builtin_whittaker_m(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    Expr** args = res->data.function.args;
    if (argc != 3)
        return builtin_arg_error("WhittakerM", argc, 3, 3);
    Expr* nd = whit_ndarray(sf_machine_whittaker_m, sf_machine_whittaker_m_c,
                            args[0], args[1], args[2], res);
    if (nd) return nd;
    Expr* m = whit_machine(res, sf_machine_whittaker_m, sf_machine_whittaker_m_c);
    if (m) return m;
    return whittaker_m_eval(args[0], args[1], args[2]);
}

/* ------------------------------------------------------------------ */
/* Registration                                                        */
/* ------------------------------------------------------------------ */

void whittaker_init(void) {
    symtab_add_builtin("WhittakerM", builtin_whittaker_m);
    symtab_get_def("WhittakerM")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    symtab_add_builtin("WhittakerW", builtin_whittaker_w);
    symtab_get_def("WhittakerW")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstrings live in info.c (info_init). */
}
