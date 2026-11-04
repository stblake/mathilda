/* Mathilda -- MathieuC[a, q, z] and MathieuCPrime[a, q, z].
 *
 *   MathieuC[a, q, z]       the even Mathieu function with characteristic value
 *                           a and parameter q: the even solution of
 *                           y'' + (a - 2 q cos 2z) y = 0.
 *   MathieuCPrime[a, q, z]  its z-derivative.
 *
 * Definition and normalization (sf_machine_mathieu.h): the Floquet form
 * sum_k c_k cos((nu + 2k) z), nu the characteristic exponent, with
 *   generic nu   sum_k c_k^2 = 1 (analytic in a and q), sign Re(c_0) > 0;
 *   integer nu   (a a characteristic value) the DLMF normalization;
 *   q == 0       cos(sqrt(a) z), and 1/sqrt 2 for a == 0.
 * For real a, q in a stability band (nu real) and for every characteristic value
 * this coincides with the reference CAS; in the instability bands and for complex
 * a or q (complex nu) the reference uses a non-analytic normalization, and the
 * two differ by a z-independent constant factor (both are the even solution).
 *
 * Evaluation, cheapest first:
 *
 *   exact rules (any arguments)
 *     MathieuC[a, 0, z]        -> Cos[Sqrt[a] z]    (1/Sqrt[2] for a == 0)
 *     MathieuCPrime[a, 0, z]   -> -Sqrt[a] Sin[Sqrt[a] z]   (0 for a == 0)
 *     MathieuCPrime[a, q, 0]   -> 0                  (ce is even)
 *     MathieuC[a, q, -z]       -> MathieuC[a, q, z]  (even, superficial sign)
 *     MathieuCPrime[a, q, -z]  -> -MathieuCPrime[a, q, z]
 *   a visible NDArray / packed buffer z   -> the machine helper threaded over it
 *   machine-precision call                -> sf_machine_mathieu (the SAME helper
 *                                            the buffer loop, the Compile VM and
 *                                            the NDArray kernel call: bit-for-bit)
 *   machine declined / arbitrary precision -> flint_num_mathieu (arb ball kernel
 *                                            in a precision ladder); a machine
 *                                            call is demoted back to machine
 *   all exact                             -> unevaluated
 *
 * Attributes: Listable, NHoldFirst, NumericFunction, Protected. Docstrings live in
 * info.c. The D[] rule lives in calculus/deriv.c; Series follows from it.
 */
#include "mathieu.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#include "arithmetic.h"          /* is_complex, make_complex, expr_is_superficially_negative */
#include "attr.h"
#include "common.h"              /* builtin_arg_error, common_machine_real_value */
#include "eval.h"
#include "expr.h"
#include "flint_num_bridge.h"    /* flint_num_mathieu */
#include "ndarray.h"
#include "numeric.h"             /* numeric_min_inexact_bits */
#include "sf_machine_mathieu.h"
#include "sf_scalar.h"           /* sf_machine_scalar_args, sf_machine_leaf */
#include "symtab.h"
#include "../sym_names.h"

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* ------------------------------------------------------------------ */
/* Small builders                                                      */
/* ------------------------------------------------------------------ */

static Expr* mk_fn1(const char* h, Expr* a) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* h, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b }, 2);
}
static Expr* mk_fn3_sym(const char* h, Expr* a, Expr* b, Expr* c) {
    return expr_new_function(expr_new_symbol(h), (Expr*[]){ a, b, c }, 3);
}
static Expr* eval_owned(Expr* e) {
    Expr* r = evaluate(e);
    expr_free(e);
    return r;
}
static bool is_exact_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Minimum precision (bits) of e's inexact parts; 0 when exact. */
static long mth_inexact_bits(const Expr* e) {
#ifdef USE_MPFR
    return numeric_min_inexact_bits(e);
#else
    Expr *re, *im;
    if (e && e->type == EXPR_REAL) return 53;
    if (e && is_complex((Expr*)e, &re, &im)
            && (re->type == EXPR_REAL || im->type == EXPR_REAL)) return 53;
    return 0;
#endif
}

/* A real numeric leaf (machine / MPFR Real, Integer, Rational) as a double. */
static bool leaf_d(const Expr* e, double* out) {
#ifdef USE_MPFR
    if (e && e->type == EXPR_MPFR) { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    return common_machine_real_value(e, out);
}

/* Demote an arb result computed at 53 bits to machine numbers. Owns `e`. */
static Expr* demote_machine(Expr* e) {
    double d;
    if (leaf_d(e, &d)) { expr_free(e); return expr_new_real(d); }
    Expr *re, *im;
    double dr, di;
    if (is_complex(e, &re, &im) && leaf_d(re, &dr) && leaf_d(im, &di)) {
        expr_free(e);
        return make_complex(expr_new_real(dr), expr_new_real(di));
    }
    return e;
}

/* ------------------------------------------------------------------ */
/* Machine scalar and buffer paths                                     */
/* ------------------------------------------------------------------ */

/* Build the machine result: Real for a real value from real inputs, else
 * Complex (a complex input gives a complex output). */
static Expr* machine_result(double rr, double ri, bool any_complex_in) {
    if (!any_complex_in && ri == 0.0) return expr_new_real(rr);
    return sf_machine_complex_result(rr, ri);
}

static Expr* mathieu_machine(Expr* res, int deriv, bool* declined) {
    SFMachineArgs ma;
    *declined = false;
    if (!sf_machine_scalar_args(res, &ma) || ma.n != 3) return NULL;
    double rr, ri;
    if (!sf_machine_mathieu(deriv, ma.re[0], ma.im[0], ma.re[1], ma.im[1],
                            ma.re[2], ma.im[2], &rr, &ri)) {
        *declined = true;
        return NULL;
    }
    return machine_result(rr, ri, ma.any_complex);
}

/* The buffer path: z a visible NDArray / packed List, a and q constant machine
 * leaves (or exact). Anything else -- a buffer a or q, an int64 / bool buffer, a
 * declined element -- goes through ndarray_delist_and_reeval, where each element
 * takes the scalar route (never NULL once a buffer is present: the scalar driver
 * cannot consume one). */
static Expr* mathieu_ndarray(Expr* res, int deriv) {
    Expr** args = res->data.function.args;
    const Expr* a_e = args[0];
    const Expr* q_e = args[1];
    const Expr* z_e = args[2];
    if (is_ndarray(a_e) || is_ndarray(q_e)) return ndarray_delist_and_reeval(res);
    if (!is_ndarray(z_e)) return NULL;
    NDType dt = z_e->data.ndarray.dtype;
    double ad, qd; bool ixa = false, ixq = false;
    if (dt == NDT_INT64 || dt == NDT_BOOL ||
        !sf_machine_leaf(a_e, &ad, &ixa) || !sf_machine_leaf(q_e, &qd, &ixq))
        return ndarray_delist_and_reeval(res);
    bool zc = ndt_is_complex(dt);
    size_t sz = ndarray_size(z_e);
    if (dt == NDT_FLOAT64) {                    /* hot path: real in, real out */
        const double* x = (const double*)z_e->data.ndarray.data;
        double* o = malloc(sizeof(double) * (sz ? sz : 1));
        if (!o) return ndarray_delist_and_reeval(res);
        size_t i = 0;
        for (; i < sz; i++) {
            double ri;
            if (!sf_machine_mathieu(deriv, ad, 0.0, qd, 0.0, x[i], 0.0, &o[i], &ri) ||
                ri != 0.0) break;
        }
        if (i == sz)
            return expr_new_ndarray_like(z_e, z_e->data.ndarray.rank,
                                         z_e->data.ndarray.dims, o, NDT_FLOAT64);
        free(o);                                /* complex-valued or declined: general */
    }
    double* buf = malloc(sizeof(double) * 2 * (sz ? sz : 1));
    if (!buf) return ndarray_delist_and_reeval(res);
    bool any_c = zc;
    for (size_t i = 0; i < sz; i++) {
        double xr, xi, rr, ri;
        ndt_get(z_e->data.ndarray.data, i, dt, &xr, &xi);
        if (!sf_machine_mathieu(deriv, ad, 0.0, qd, 0.0, xr, xi, &rr, &ri)) {
            free(buf);
            return ndarray_delist_and_reeval(res);
        }
        buf[2 * i] = rr; buf[2 * i + 1] = ri;
        if (ri != 0.0) any_c = true;
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
/* Exact rules                                                         */
/* ------------------------------------------------------------------ */

static Expr* mathieu_exact(const char* head, Expr* a, Expr* q, Expr* z, int deriv,
                           bool any_inexact) {
    /* q == 0: the Mathieu equation is y'' + a y = 0. */
    if (is_exact_zero(q)) {
        /* the constant answers (1/Sqrt[2], 0) are exact, so they are only
         * returned for exact input; an inexact z takes the numeric path */
        if (is_exact_zero(a) && any_inexact) return NULL;
        if (is_exact_zero(a))
            return deriv ? expr_new_integer(0)
                         : eval_owned(mk_fn2("Power", expr_new_integer(2),
                                             mk_fn2("Rational", expr_new_integer(-1),
                                                    expr_new_integer(2))));
        Expr* sa = mk_fn1("Sqrt", expr_copy(a));
        Expr* arg = mk_fn2("Times", expr_copy(sa), expr_copy(z));
        if (!deriv) { expr_free(sa); return eval_owned(mk_fn1("Cos", arg)); }
        return eval_owned(expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            expr_new_integer(-1), sa, mk_fn1("Sin", arg) }, 3));
    }
    /* ce'(0) = 0 (ce even); an inexact a or q gives the numeric 0. */
    if (deriv && is_exact_zero(z) && !any_inexact) return expr_new_integer(0);
    /* Parity in z, on the superficial sign (as Cos / Sin do). */
    if (expr_is_superficially_negative(z)) {
        Expr* nz = eval_owned(mk_fn2("Times", expr_new_integer(-1), expr_copy(z)));
        Expr* f = mk_fn3_sym(head, expr_copy(a), expr_copy(q), nz);
        if (!deriv) return f;
        return mk_fn2("Times", expr_new_integer(-1), f);
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
/* Driver                                                              */
/* ------------------------------------------------------------------ */

static Expr* mathieu_builtin(Expr* res, const char* head, int deriv) {
    if (res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc != 3) return builtin_arg_error(head, argc, 3, 3);
    Expr** args = res->data.function.args;

    Expr* nd = mathieu_ndarray(res, deriv);
    if (nd) return nd;

    long b0 = mth_inexact_bits(args[0]), b1 = mth_inexact_bits(args[1]);
    long b2 = mth_inexact_bits(args[2]);
    Expr* ex = mathieu_exact(head, args[0], args[1], args[2], deriv,
                             b0 > 0 || b1 > 0 || b2 > 0);
    if (ex) return ex;

    if (b0 == 0 && b1 == 0 && b2 == 0) return NULL;          /* all exact: stay */
    if (!expr_is_numeric_quantity(args[0]) || !expr_is_numeric_quantity(args[1]) ||
        !expr_is_numeric_quantity(args[2])) return NULL;      /* symbolic */
    long outb = 0;
    if (b0 > 0) outb = b0;
    if (b1 > 0 && (outb == 0 || b1 < outb)) outb = b1;
    if (b2 > 0 && (outb == 0 || b2 < outb)) outb = b2;

    if (outb <= 53) {
        bool declined = false;
        Expr* m = mathieu_machine(res, deriv, &declined);
        if (m) return m;
        /* Machine declined (cancellation at large Im z, a large |a| or |q|, a
         * degenerate normalization): the ball kernel at 53 bits, demoted. */
        Expr* v = flint_num_mathieu(args[0], args[1], args[2], deriv);
        if (v) return demote_machine(v);
        return NULL;
    }
    return flint_num_mathieu(args[0], args[1], args[2], deriv);
}

Expr* builtin_mathieuc(Expr* res)      { return mathieu_builtin(res, "MathieuC", 0); }
Expr* builtin_mathieucprime(Expr* res) { return mathieu_builtin(res, "MathieuCPrime", 1); }

void mathieu_init(void) {
    const unsigned attrs = ATTR_LISTABLE | ATTR_NHOLDFIRST | ATTR_NUMERICFUNCTION |
                           ATTR_PROTECTED;
    symtab_add_builtin("MathieuC", builtin_mathieuc);
    symtab_get_def("MathieuC")->attributes |= attrs;
    symtab_add_builtin("MathieuCPrime", builtin_mathieucprime);
    symtab_get_def("MathieuCPrime")->attributes |= attrs;
    /* Docstrings live in info.c (info_init). */
}
