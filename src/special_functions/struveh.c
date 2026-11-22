/* Mathilda -- the Struve function H_nu(z).
 *
 *   StruveH[nu, z]   H_nu(z) = (z/2)^(nu+1) (2/(Sqrt[Pi] Gamma[nu+3/2]))
 *                              1F2(1; 3/2, nu+3/2; -z^2/4)
 *
 * Evaluation:
 *   StruveH[nu, 0] = 0            exact, for a real nu > -1
 *   machine precision, real nu,z  the shared double kernel sf_machine_struve_h
 *                                 (bit-for-bit with NDArray / Compile[])
 *   anything else numeric         arb's regularized 1F2 (flint_num_struveh),
 *                                 demoted to machine numbers for machine input;
 *                                 the unregularized 1F2 template only without FLINT
 *   exact half-odd-integer nu     elementary closed form (sh_half_integer)
 *   exact integer nu = -n < 0     (-1)^n H_n(z) + a Laurent polynomial
 *                                 (struve_negative_integer_order)
 *   everything else               stays symbolic (return NULL)
 *
 * The 1F2 identity is also the reduction used by the half-line Mellin transform
 * (reduce_to_hypergeometric in integrate_ramanujan.c), so Integrate[x^(s-1)
 * StruveH[0, a x], {x,0,Infinity}] closes through rec_pfq.
 *
 * Attributes: Listable, NumericFunction, Protected.
 */
#include "struveh.h"

#include "symtab.h"
#include "attr.h"
#include "eval.h"
#include "parse.h"
#include "sym_names.h"
#include "common.h"
#include "numeric.h"
#include "arithmetic.h"         /* is_complex */
#include "flint_num_bridge.h"   /* flint_num_struveh (arb, arbitrary precision) */
#include "sf_machine.h"         /* sf_machine_struve_h (shared double kernel)   */

#include <stdbool.h>
#include <stdlib.h>

#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* Small expr builders (local idiom). */
static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_fn1(const char* name, Expr* a) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a }, 1);
}
static Expr* mk_fn2(const char* name, Expr* a, Expr* b) {
    return expr_new_function(expr_new_symbol(name), (Expr*[]){ a, b }, 2);
}
/* evaluate `call` to a fixed point, then free it (evaluate borrows its arg). */
static Expr* sh_eval(Expr* call) {
    Expr* r = evaluate(call);
    expr_free(call);
    return r;
}

/* H_nu(z) via its 1F2 representation, numericised; NULL if the result is not a
 * number (e.g. the pFq evaluator declined). */
static Expr* struveh_numeric(const Expr* nu, const Expr* z) {
    static const char* TMPL =
        "(zz/2)^(nn+1) (2/(Sqrt[Pi] Gamma[nn+3/2])) "
        "HypergeometricPFQ[{1}, {3/2, nn+3/2}, -zz^2/4]";
    Expr* tmpl = parse_expression(TMPL);
    if (!tmpl) return NULL;
    Expr* rules = expr_new_function(expr_new_symbol("List"), (Expr*[]){
        mk_fn2("Rule", expr_new_symbol("nn"), expr_copy((Expr*)nu)),
        mk_fn2("Rule", expr_new_symbol("zz"), expr_copy((Expr*)z)) }, 2);
    Expr* subbed = sh_eval(mk_fn2("ReplaceAll", tmpl, rules));
    Expr* r = sh_eval(mk_fn1("N", subbed));
    if (r && expr_is_numeric_quantity(r)) return r;
    if (r) expr_free(r);
    return NULL;
}

/* A machine Real, or Complex of machine Reals. */
static bool sh_is_machine_number(const Expr* e) {
    Expr *re, *im;
    if (!e) return false;
    if (e->type == EXPR_REAL) return true;
    return is_complex((Expr*)e, &re, &im) &&
           re->type == EXPR_REAL && im->type == EXPR_REAL;
}

/* A real numeric leaf (machine/MPFR Real, Integer, Rational) as a double. */
static bool sh_leaf_d(const Expr* e, double* out) {
#ifdef USE_MPFR
    if (e && e->type == EXPR_MPFR) { *out = mpfr_get_d(e->data.mpfr, MPFR_RNDN); return true; }
#endif
    return common_machine_real_value(e, out);
}

/* Demote an arb result (53-bit MPFR leaves, possibly Complex) to machine Reals.
 * Takes ownership of `e`; returns it unchanged if it is not a plain number. */
static Expr* sh_demote_machine(Expr* e) {
    double d;
    if (sh_leaf_d(e, &d)) { expr_free(e); return expr_new_real(d); }
    Expr *re, *im;
    double dr, di;
    if (is_complex(e, &re, &im) && sh_leaf_d(re, &dr) && sh_leaf_d(im, &di)) {
        expr_free(e);
        return make_complex(expr_new_real(dr), expr_new_real(di));
    }
    return e;
}

/* Minimum precision (bits) of e's inexact components (machine Real / MPFR,
 * or Complex of those); 0 when e is exact, i.e. N[] has not been applied. */
static long sh_inexact_bits(const Expr* e) {
#ifdef USE_MPFR
    return numeric_min_inexact_bits(e);
#else
    /* numeric_min_inexact_bits is MPFR-only; without MPFR the only inexact
     * leaf is a machine Real (53 bits), possibly inside Complex[...]. */
    Expr *re, *im;
    if (e && e->type == EXPR_REAL) return 53;
    if (e && is_complex((Expr*)e, &re, &im)
            && (re->type == EXPR_REAL || im->type == EXPR_REAL)) return 53;
    return 0;
#endif
}
static bool sh_inexact(const Expr* e) { return sh_inexact_bits(e) > 0; }

/* ---- half-odd-integer order: elementary closed forms ------------------- */

static Expr* sh_rat(int64_t p, int64_t q) {
    return mk_fn2("Rational", mk_int(p), mk_int(q));
}

/* a_k(n + 1/2) = (n + k)! / (2^k k! (n - k)!), the Hankel coefficient
 * (DLMF 10.49.1), as an unevaluated exact expression. */
static Expr* sh_hankel_a(int64_t n, int64_t k) {
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        mk_fn1("Factorial", mk_int(n + k)),
        mk_fn2("Power", mk_int(2), mk_int(-k)),
        mk_fn2("Power", mk_fn1("Factorial", mk_int(k)), mk_int(-1)),
        mk_fn2("Power", mk_fn1("Factorial", mk_int(n - k)), mk_int(-1)) }, 4);
}

/* Sum_{k = par, par + 2, ... <= n} (-1)^((k - par)/2) a_k(n + 1/2) z^-k. */
static Expr* sh_hankel_sum(int64_t n, int par, const Expr* z) {
    int64_t cnt = n >= par ? (n - par) / 2 + 1 : 0;
    if (cnt == 0) return mk_int(0);
    Expr** t = malloc(sizeof(Expr*) * (size_t)cnt);
    if (!t) return NULL;
    for (int64_t i = 0; i < cnt; i++) {
        int64_t k = par + 2 * i;
        t[i] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_int((i & 1) ? -1 : 1), sh_hankel_a(n, k),
            mk_fn2("Power", expr_copy((Expr*)z), mk_int(-k)) }, 3);
    }
    Expr* s = expr_new_function(expr_new_symbol("Plus"), t, (size_t)cnt);
    free(t);
    return s;
}

/* The spherical-Bessel form of J_{n+1/2}(z) (want_y false) or Y_{n+1/2}(z)
 * (DLMF 10.49.2/10.49.4, times Sqrt(2 z / Pi)):
 *
 *   J = Sqrt(2/Pi) z^(-1/2) ( S E(z) + C O(z))
 *   Y = Sqrt(2/Pi) z^(-1/2) (-C E(z) + S O(z))
 *
 * with S = Sin(z - n Pi/2), C = Cos(z - n Pi/2) (reduced to +-Sin z, +-Cos z),
 * E the even and O the odd Hankel sum. */
static Expr* sh_sph_bessel(int64_t n, bool want_y, const Expr* z) {
    Expr* E = sh_hankel_sum(n, 0, z);
    Expr* O = sh_hankel_sum(n, 1, z);
    if (!E || !O) { expr_free(E); expr_free(O); return NULL; }
    int r = (int)(n % 4);
    /* sin(z - r Pi/2), cos(z - r Pi/2) */
    static const int s_sgn[4] = { 1, -1, -1, 1 }, s_sin[4] = { 1, 0, 1, 0 };
    static const int c_sgn[4] = { 1, 1, -1, -1 };
    Expr* S = mk_fn2("Times", mk_int(s_sgn[r]),
                     mk_fn1(s_sin[r] ? "Sin" : "Cos", expr_copy((Expr*)z)));
    Expr* C = mk_fn2("Times", mk_int(c_sgn[r]),
                     mk_fn1(s_sin[r] ? "Cos" : "Sin", expr_copy((Expr*)z)));
    Expr* body = want_y
        ? mk_fn2("Plus", expr_new_function(expr_new_symbol("Times"),
                     (Expr*[]){ mk_int(-1), C, E }, 3), mk_fn2("Times", S, O))
        : mk_fn2("Plus", mk_fn2("Times", S, E), mk_fn2("Times", C, O));
    return expr_new_function(expr_new_symbol("Times"), (Expr*[]){
        mk_fn2("Power", mk_int(2), sh_rat(1, 2)),
        mk_fn2("Power", expr_new_symbol("Pi"), sh_rat(-1, 2)),
        mk_fn2("Power", expr_copy((Expr*)z), sh_rat(-1, 2)),
        body }, 4);
}

/* A half-odd-integer order nu = k/2 (k odd): the elementary closed forms
 * (DLMF 11.4.4 / 11.4.5), with n >= 0:
 *
 *   H_{-(n+1/2)}(z) = (-1)^n J_{n+1/2}(z)
 *   H_{n+1/2}(z)    = Y_{n+1/2}(z)
 *                     + (1/Pi) Sum_{j=0}^n Gamma(j+1/2) (z/2)^(n-2j-1/2) / (n-j)!
 *
 * with the half-integer Bessel functions written out in spherical form, so the
 * result is elementary for any n (no Bessel DownValues are involved). */
static Expr* sh_half_integer(const Expr* nu, const Expr* z) {
    if (nu->type != EXPR_FUNCTION || !is_rational_like(nu)) return NULL;
    const Expr* num = nu->data.function.args[0];
    const Expr* den = nu->data.function.args[1];
    if (num->type != EXPR_INTEGER || den->type != EXPR_INTEGER || den->data.integer != 2)
        return NULL;
    int64_t k = num->data.integer;                 /* odd: Rational is canonical */
    int64_t n = k < 0 ? (-k - 1) / 2 : (k - 1) / 2;
    Expr* bes = sh_sph_bessel(n, k > 0, z);
    if (!bes) return NULL;
    if (k < 0)
        return sh_eval(mk_fn2("Times", mk_int((n & 1) ? -1 : 1), bes));
    Expr** t = malloc(sizeof(Expr*) * (size_t)(n + 1));
    if (!t) { expr_free(bes); return NULL; }
    for (int64_t j = 0; j <= n; j++) {
        /* Gamma(j+1/2) 2^-(n-2j-1/2) z^(n-2j-1/2) / (n-j)! */
        Expr* e = sh_rat(2 * (n - 2 * j) - 1, 2);
        t[j] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_fn1("Gamma", sh_rat(2 * j + 1, 2)),
            mk_fn2("Power", mk_int(2), mk_fn2("Times", mk_int(-1), expr_copy(e))),
            mk_fn2("Power", expr_copy((Expr*)z), e),
            mk_fn2("Power", mk_fn1("Factorial", mk_int(n - j)), mk_int(-1)) }, 4);
    }
    Expr* poly = expr_new_function(expr_new_symbol("Plus"), t, (size_t)(n + 1));
    free(t);
    Expr* part = mk_fn2("Times", mk_fn2("Power", expr_new_symbol("Pi"), mk_int(-1)), poly);
    return sh_eval(mk_fn2("Plus", part, bes));
}

/* ---- negative integer order: reflection to the order +n ---------------
 *
 * The defining series (DLMF 11.2.1 / 11.2.2)
 *
 *   H_nu(z) = Sum_{k>=0} (-1)^k (z/2)^(2k+nu+1) / (Gamma(k+3/2) Gamma(k+nu+3/2))
 *   L_nu(z) = Sum_{k>=0}        (z/2)^(2k+nu+1) / (Gamma(k+3/2) Gamma(k+nu+3/2))
 *
 * has no Gamma pole at nu = -n (integer n >= 1): k - n + 3/2 is a half-odd
 * integer.  Re-indexing its terms k >= n by k = j + n gives term-for-term
 * (-1)^n times (H) or exactly (L) the series of order +n, so
 *
 *   H_{-n}(z) = (-1)^n H_n(z) + Sum_{k=0}^{n-1} (-1)^k (z/2)^(2k-n+1) / (Gamma(k+3/2) Gamma(k-n+3/2))
 *   L_{-n}(z) =        L_n(z) + Sum_{k=0}^{n-1}        (z/2)^(2k-n+1) / (Gamma(k+3/2) Gamma(k-n+3/2))
 *
 * e.g. H_{-1} = 2/Pi - H_1 and L_{-1} = 2/Pi + L_1 (consistent with the
 * recurrence H_{nu-1} + H_{nu+1} = (2 nu/z) H_nu + (z/2)^nu/(Sqrt(Pi) Gamma(nu+3/2))
 * at nu = 0).  The finite sum is a Laurent polynomial with exact coefficients.
 * Returns NULL unless nu is an exact Integer <= -1 (and n below a sanity cap,
 * since the sum has n terms). */
Expr* struve_negative_integer_order(const Expr* nu, const Expr* z, bool modified) {
    if (!nu || nu->type != EXPR_INTEGER || nu->data.integer >= 0) return NULL;
    if (nu->data.integer < -100000) return NULL;
    int64_t n = -nu->data.integer;
    Expr** t = malloc(sizeof(Expr*) * (size_t)(n + 1));
    if (!t) return NULL;
    const char* head = modified ? "StruveL" : "StruveH";
    t[0] = mk_fn2("Times", mk_int((!modified && (n & 1)) ? -1 : 1),
                  mk_fn2(head, mk_int(n), expr_copy((Expr*)z)));
    for (int64_t k = 0; k < n; k++) {
        /* (+-1) (z/2)^(2k-n+1) / (Gamma(k+3/2) Gamma(k-n+3/2)) */
        int64_t e = 2 * k - n + 1;
        t[k + 1] = expr_new_function(expr_new_symbol("Times"), (Expr*[]){
            mk_int((!modified && (k & 1)) ? -1 : 1),
            /* (z/2)^0 is spelled 1, so H_{-1}(0) = 2/Pi rather than 0^0. */
            e == 0 ? mk_int(1)
                   : mk_fn2("Power", mk_fn2("Times", sh_rat(1, 2), expr_copy((Expr*)z)),
                            mk_int(e)),
            mk_fn2("Power", mk_fn1("Gamma", sh_rat(2 * k + 3, 2)), mk_int(-1)),
            mk_fn2("Power", mk_fn1("Gamma", sh_rat(2 * (k - n) + 3, 2)), mk_int(-1)) }, 4);
    }
    Expr* sum = expr_new_function(expr_new_symbol("Plus"), t, (size_t)(n + 1));
    free(t);
    return sh_eval(sum);
}

Expr* builtin_struveh(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2) return NULL;
    Expr* nu = res->data.function.args[0];
    Expr* z  = res->data.function.args[1];

    /* H_nu(0) = 0 for a real nu > -1 (exact). */
    if (z->type == EXPR_INTEGER && z->data.integer == 0) {
        double nuv;
        if (common_machine_real_value(nu, &nuv) && nuv > -1.0) return mk_int(0);
    }

    /* Numeric path: at least one inexact argument, both concrete numbers. */
    if ((sh_inexact(nu) || sh_inexact(z)) &&
        expr_is_numeric_quantity(nu) && expr_is_numeric_quantity(z)) {
        /* Target precision: the min inexact-bit count over the arguments
         * (an exact arg imposes no constraint), floored implicitly at 53. */
        long bnu = sh_inexact_bits(nu);
        long bz  = sh_inexact_bits(z);
        long outb = bnu;
        if (bz > 0 && (outb == 0 || bz < outb)) outb = bz;

        if (outb > 53) {
            /* Arbitrary precision: the rigorous arb kernel (all nu, all C).
             * On a genuine decline (non-finite) leave it symbolic rather than
             * silently downgrade to machine precision. */
            Expr* v = flint_num_struveh(nu, z);
            if (v) return v;
        } else {
            /* Machine precision: the shared double kernel for real nu, z.
             * Bit-identical to the NDArray / Compile surfaces. */
            double nuv, zv;
            if (common_machine_real_value(nu, &nuv) &&
                common_machine_real_value(z, &zv)) {
                double h;
                if (sf_machine_struve_h(nuv, zv, &h)) return expr_new_real(h);
            }
            /* Kernel declined (large nu, a real z < 0 with a complex value, a
             * complex z, a non-finite double): the rigorous arb kernel at 53 bits,
             * demoted to machine numbers (machine in -> machine out). arb works
             * from the REGULARIZED 1F2, so the Gamma(nu+3/2) poles cost nothing. */
            Expr* v = flint_num_struveh(nu, z);
            if (v) return sh_demote_machine(v);
            /* No FLINT: the 1F2 template, kept only when it collapsed to a
             * machine number (it is unregularized: at a Gamma(nu+3/2) pole it
             * yields 0. * HypergeometricPFQ[...], which must not be returned). */
            if (z->type == EXPR_REAL) {
                double nuv2;
                if (common_machine_real_value(nu, &nuv2)) {
                    Expr* t = struveh_numeric(nu, z);
                    if (t && sh_is_machine_number(t)) return t;
                    expr_free(t);
                }
            }
        }
        return NULL;
    }
    /* Exact negative integer order: reflect to the order +n. */
    Expr* neg = struve_negative_integer_order(nu, z, false);
    if (neg) return neg;
    /* Exact half-odd-integer order: elementary. */
    return sh_half_integer(nu, z);
}

void struveh_init(void) {
    symtab_add_builtin("StruveH", builtin_struveh);
    symtab_get_def("StruveH")->attributes |=
        (ATTR_LISTABLE | ATTR_NUMERICFUNCTION | ATTR_PROTECTED);
    /* Docstring lives in info.c (info_init). */
}
