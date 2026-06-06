#include "complex.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include "numeric.h"
#include "common.h"
#include "sym_names.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#ifdef USE_MPFR
#include <mpfr.h>
#endif

/* True iff `e` reduces to a real machine number under N[]. Used to gate
 * complex decomposition so the structural Re/Im split is only applied
 * when the resulting real and imaginary parts are concretely numeric
 * (Mathematica's behaviour: symbolic args like `Abs[1 + I x]` stay
 * unevaluated rather than guess the reality of x). */
static bool is_numeric_real(Expr* e) {
    if (!e) return false;
    int64_t n, d;
    if (e->type == EXPR_INTEGER || e->type == EXPR_REAL) return true;
    if (is_rational(e, &n, &d)) return true;
    Expr* approx = numericalize(e, numeric_machine_spec());
    bool ok = false;
    if (approx) {
        if (approx->type == EXPR_INTEGER || approx->type == EXPR_REAL) ok = true;
        else if (is_rational(approx, &n, &d)) ok = true;
        expr_free(approx);
    }
    return ok;
}

/* Structurally decompose an expression into real and imaginary parts.
 *
 * Returns true when at least one Complex[a, b] literal is found nested
 * inside Plus or Times factors. The returned re_out/im_out are freshly
 * evaluated expressions owned by the caller.
 *
 * Subexpressions that are not Complex literals (symbols, Sqrt[..], Sin[..],
 * arbitrary user functions, etc.) are assumed to be real-valued. This is
 * the same assumption Mathematica makes when ComplexExpand is not in play
 * and matches Mathilda's existing handling of `1 + I` (which already
 * folds into Complex[1, 1]); the gap was purely the Plus/Times of Complex
 * literals with non-numeric (symbolic) coefficients.
 *
 * Relies on the evaluator's existing folding so Times has at most one
 * Complex factor by the time we see it.
 */
static bool complex_decompose(Expr* e, Expr** re_out, Expr** im_out) {
    if (!e) return false;

    Expr *cre, *cim;
    if (is_complex(e, &cre, &cim)) {
        *re_out = expr_copy(cre);
        *im_out = expr_copy(cim);
        return true;
    }

    if (head_is(e, SYM_Plus)) {
        size_t n = e->data.function.arg_count;
        Expr** re_terms = malloc(sizeof(Expr*) * n);
        Expr** im_terms = malloc(sizeof(Expr*) * n);
        size_t nr = 0, ni = 0;
        bool any_complex = false;
        for (size_t k = 0; k < n; k++) {
            Expr* term = e->data.function.args[k];
            Expr *tre, *tim;
            if (complex_decompose(term, &tre, &tim)) {
                any_complex = true;
                re_terms[nr++] = tre;
                im_terms[ni++] = tim;
            } else {
                re_terms[nr++] = expr_copy(term);
            }
        }
        if (!any_complex) {
            for (size_t k = 0; k < nr; k++) expr_free(re_terms[k]);
            free(re_terms);
            free(im_terms);
            return false;
        }
        Expr* re_sum;
        if (nr == 1) { re_sum = re_terms[0]; }
        else { re_sum = expr_new_function(expr_new_symbol("Plus"), re_terms, nr); }
        free(re_terms);
        Expr* im_sum;
        if (ni == 0) { im_sum = expr_new_integer(0); }
        else if (ni == 1) { im_sum = im_terms[0]; }
        else { im_sum = expr_new_function(expr_new_symbol("Plus"), im_terms, ni); }
        free(im_terms);
        *re_out = eval_and_free(re_sum);
        *im_out = eval_and_free(im_sum);
        return true;
    }

    if (head_is(e, SYM_Times)) {
        /* Walk factors left-to-right, maintaining a running complex product
         * (acc_re + acc_im*I). A factor that itself decomposes into a
         * non-trivial complex form contributes via complex multiplication
         * (ac - bd) + (ad + bc)*I; otherwise we treat the factor as real
         * and scale both parts by it. */
        Expr* acc_re = expr_new_integer(1);
        Expr* acc_im = expr_new_integer(0);
        bool any_complex = false;
        for (size_t k = 0; k < e->data.function.arg_count; k++) {
            Expr* fac = e->data.function.args[k];
            Expr* fre = NULL;
            Expr* fim = NULL;
            if (complex_decompose(fac, &fre, &fim)) {
                any_complex = true;
                Expr* ac_args[2] = { expr_copy(acc_re), expr_copy(fre) };
                Expr* ac = expr_new_function(expr_new_symbol("Times"), ac_args, 2);
                Expr* bd_args[2] = { expr_copy(acc_im), expr_copy(fim) };
                Expr* bd = expr_new_function(expr_new_symbol("Times"), bd_args, 2);
                Expr* neg_bd_args[2] = { expr_new_integer(-1), bd };
                Expr* neg_bd = expr_new_function(expr_new_symbol("Times"), neg_bd_args, 2);
                Expr* re_sum_args[2] = { ac, neg_bd };
                Expr* re_new = expr_new_function(expr_new_symbol("Plus"), re_sum_args, 2);

                Expr* ad_args[2] = { expr_copy(acc_re), expr_copy(fim) };
                Expr* ad = expr_new_function(expr_new_symbol("Times"), ad_args, 2);
                Expr* bc_args[2] = { expr_copy(acc_im), expr_copy(fre) };
                Expr* bc = expr_new_function(expr_new_symbol("Times"), bc_args, 2);
                Expr* im_sum_args[2] = { ad, bc };
                Expr* im_new = expr_new_function(expr_new_symbol("Plus"), im_sum_args, 2);

                expr_free(acc_re);
                expr_free(acc_im);
                expr_free(fre);
                expr_free(fim);
                acc_re = eval_and_free(re_new);
                acc_im = eval_and_free(im_new);
            } else {
                Expr* re_args[2] = { acc_re, expr_copy(fac) };
                Expr* re_new = expr_new_function(expr_new_symbol("Times"), re_args, 2);
                Expr* im_args[2] = { acc_im, expr_copy(fac) };
                Expr* im_new = expr_new_function(expr_new_symbol("Times"), im_args, 2);
                acc_re = eval_and_free(re_new);
                acc_im = eval_and_free(im_new);
            }
        }
        if (!any_complex) {
            expr_free(acc_re);
            expr_free(acc_im);
            return false;
        }
        *re_out = acc_re;
        *im_out = acc_im;
        return true;
    }

    return false;
}

void complex_init(void) {
    symtab_add_builtin("Re", builtin_re);
    symtab_add_builtin("Im", builtin_im);
    symtab_add_builtin("ReIm", builtin_reim);
    symtab_add_builtin("Abs", builtin_abs);
    symtab_add_builtin("Sign", builtin_sign);
    symtab_add_builtin("Conjugate", builtin_conjugate);
    symtab_add_builtin("Arg", builtin_arg);

    symtab_get_def("Re")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Im")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("ReIm")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Abs")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Sign")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Conjugate")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Arg")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);

    symtab_set_docstring("Sign",
        "Sign[x] gives -1, 0, or 1 for real numeric x according to its sign, "
        "and z/Abs[z] for a nonzero numeric complex z.");
    symtab_set_docstring("Re",
        "Re[z] gives the real part of numeric z; Re[Re[z]], Re[Im[z]], "
        "Re[Abs[z]], Re[Arg[z]] fold since those heads are real-valued.");
    symtab_set_docstring("Im",
        "Im[z] gives the imaginary part of numeric z, and 0 for real or "
        "real-valued (Re/Im/Abs/Arg) arguments.");
    symtab_set_docstring("Abs",
        "Abs[z] gives the absolute value (modulus) of numeric z, "
        "Sqrt[Re[z]^2 + Im[z]^2] for complex z.");
    symtab_set_docstring("Arg",
        "Arg[z] gives the argument (phase angle in (-Pi, Pi]) of numeric z; "
        "0 for nonnegative reals, Pi for negative reals.");
    symtab_set_docstring("Conjugate",
        "Conjugate[z] gives the complex conjugate Re[z] - I Im[z] of numeric z; "
        "real and real-valued (Re/Im/Abs/Arg) arguments are returned unchanged.");
    symtab_set_docstring("ReIm",
        "ReIm[z] gives {Re[z], Im[z]}, the real and imaginary parts of numeric z "
        "as a list; real-valued arguments give {z, 0}.");
}

/* True when `e` is `f[x]` for one of the real-valued-by-construction heads
 * Re, Im, Abs, Arg. Used by Re/Im to fold Re[Re[z]] -> Re[z], Im[Abs[z]] -> 0,
 * etc. Mirrors the same set Conjugate treats as fixed points. */
static bool is_real_valued_head_call(Expr* e) {
    return (head_is(e, SYM_Re) || head_is(e, SYM_Im) ||
            head_is(e, SYM_Abs) || head_is(e, SYM_Arg)) &&
           e->data.function.arg_count == 1;
}

Expr* builtin_re(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    /* Re[f[z]] -> f[z] when f is real-valued by construction (Re, Im, Abs, Arg). */
    if (is_real_valued_head_call(arg)) {
        res->data.function.args[0] = NULL;
        return arg;
    }
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        return expr_copy(re);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_copy(arg);
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        return expr_copy(arg);
    }
#endif
    if (complex_decompose(arg, &re, &im)) {
        if (is_numeric_real(re) && is_numeric_real(im)) {
            expr_free(im);
            return re;
        }
        expr_free(re);
        expr_free(im);
    }
    return NULL;
}

Expr* builtin_im(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    /* Im[f[z]] -> 0 when f is real-valued by construction (Re, Im, Abs, Arg). */
    if (is_real_valued_head_call(arg)) {
        return expr_new_integer(0);
    }
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        return expr_copy(im);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_new_integer(0);
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        return expr_new_integer(0);
    }
#endif
    if (complex_decompose(arg, &re, &im)) {
        if (is_numeric_real(re) && is_numeric_real(im)) {
            expr_free(re);
            return im;
        }
        expr_free(re);
        expr_free(im);
    }
    return NULL;
}

Expr* builtin_reim(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    /* ReIm[f[z]] -> {f[z], 0} when f is real-valued by construction. */
    if (is_real_valued_head_call(arg)) {
        Expr** results = malloc(sizeof(Expr*) * 2);
        res->data.function.args[0] = NULL;
        results[0] = arg;
        results[1] = expr_new_integer(0);
        Expr* list = expr_new_function(expr_new_symbol("List"), results, 2);
        free(results);
        return list;
    }
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        Expr** results = malloc(sizeof(Expr*) * 2);
        results[0] = expr_copy(re);
        results[1] = expr_copy(im);
        return expr_new_function(expr_new_symbol("List"), results, 2);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        Expr** results = malloc(sizeof(Expr*) * 2);
        results[0] = expr_copy(arg);
        results[1] = expr_new_integer(0);
        return expr_new_function(expr_new_symbol("List"), results, 2);
    }
    if (complex_decompose(arg, &re, &im)) {
        if (is_numeric_real(re) && is_numeric_real(im)) {
            Expr** results = malloc(sizeof(Expr*) * 2);
            results[0] = re;
            results[1] = im;
            Expr* list = expr_new_function(expr_new_symbol("List"), results, 2);
            free(results);
            return list;
        }
        expr_free(re);
        expr_free(im);
    }
    return NULL;
}

Expr* builtin_abs(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    bool from_complex = is_complex(arg, &re, &im);
    Expr* re_owned = NULL;
    Expr* im_owned = NULL;
    if (from_complex) {
        re_owned = expr_copy(re);
        im_owned = expr_copy(im);
    } else if (complex_decompose(arg, &re_owned, &im_owned)) {
        if (is_numeric_real(re_owned) && is_numeric_real(im_owned)) {
            from_complex = true;
        } else {
            expr_free(re_owned);
            expr_free(im_owned);
            re_owned = NULL;
            im_owned = NULL;
        }
    }
    if (from_complex) {
#ifdef USE_MPFR
        /* MPFR fast path: when at least one component carries MPFR, fold
         * directly through mpfr_hypot rather than constructing the
         * symbolic Sqrt[Plus[Power[re,2], Power[im,2]]] tree. The
         * symbolic form would otherwise re-enter the MPFR Power /
         * Plus chain, allocating two squarings and a sum at each
         * intermediate step. mpfr_hypot is also stable when |re| and
         * |im| span very different magnitudes. */
        if (numeric_any_mpfr(re_owned, im_owned)) {
            long bits = numeric_combined_bits(re_owned, im_owned, 0);
            mpfr_t a_re, a_im, b_re, b_im;
            mpfr_init2(a_re, bits); mpfr_init2(a_im, bits);
            mpfr_init2(b_re, bits); mpfr_init2(b_im, bits);
            bool ok_r = get_approx_mpfr(re_owned, a_re, a_im, NULL);
            bool ok_i = get_approx_mpfr(im_owned, b_re, b_im, NULL);
            if (ok_r && ok_i && mpfr_zero_p(a_im) && mpfr_zero_p(b_im)) {
                Expr* result = expr_new_mpfr_bits(bits);
                mpfr_hypot(result->data.mpfr, a_re, b_re, MPFR_RNDN);
                mpfr_clear(a_re); mpfr_clear(a_im);
                mpfr_clear(b_re); mpfr_clear(b_im);
                expr_free(re_owned);
                expr_free(im_owned);
                return result;
            }
            mpfr_clear(a_re); mpfr_clear(a_im);
            mpfr_clear(b_re); mpfr_clear(b_im);
        }
#endif
        Expr* pow_re_args[2] = { re_owned, expr_new_integer(2) };
        Expr* pow_re = expr_new_function(expr_new_symbol("Power"), pow_re_args, 2);

        Expr* pow_im_args[2] = { im_owned, expr_new_integer(2) };
        Expr* pow_im = expr_new_function(expr_new_symbol("Power"), pow_im_args, 2);

        Expr* sum_args[2] = { pow_re, pow_im };
        Expr* sum = expr_new_function(expr_new_symbol("Plus"), sum_args, 2);

        Expr* rat_args[2] = { expr_new_integer(1), expr_new_integer(2) };
        Expr* half = expr_new_function(expr_new_symbol("Rational"), rat_args, 2);
        Expr* pow_args[2] = { sum, half };
        return expr_new_function(expr_new_symbol("Power"), pow_args, 2);
    }
    int64_t n, d;
    if (arg->type == EXPR_BIGINT) {
        mpz_t r;
        mpz_init(r);
        mpz_abs(r, arg->data.bigint);
        Expr* result = expr_bigint_normalize(expr_new_bigint_from_mpz(r));
        mpz_clear(r);
        return result;
    }
    if (arg->type == EXPR_INTEGER) {
        int64_t val = arg->data.integer;
        return expr_new_integer(val < 0 ? -val : val);
    }
    if (arg->type == EXPR_REAL) {
        double val = arg->data.real;
        return expr_new_real(val < 0.0 ? -val : val);
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        mpfr_prec_t prec = mpfr_get_prec(arg->data.mpfr);
        Expr* result = expr_new_mpfr_bits(prec);
        mpfr_abs(result->data.mpfr, arg->data.mpfr, MPFR_RNDN);
        return result;
    }
#endif
    if (is_rational(arg, &n, &d)) {
        return make_rational(n < 0 ? -n : n, d);
    }
    return NULL;
}

/* Sign[x] returns -1, 0, or 1 for a real numeric x and z/Abs[z] for a
 * nonzero numeric complex z. Symbolic or partially-numeric inputs leave
 * the call unevaluated so the symbolic head flows through the evaluator.
 *
 * Real numeric inputs reduce to an exact Integer regardless of the input
 * type (Integer/BigInt/Real/MPFR/Rational), matching Mathematica: the
 * sign is an exact quantity even when the magnitude is inexact. */
Expr* builtin_sign(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];

    if (arg->type == EXPR_INTEGER) {
        int64_t v = arg->data.integer;
        return expr_new_integer(v < 0 ? -1 : v > 0 ? 1 : 0);
    }
    if (arg->type == EXPR_BIGINT) {
        return expr_new_integer(mpz_sgn(arg->data.bigint));
    }
    if (arg->type == EXPR_REAL) {
        double v = arg->data.real;
        return expr_new_integer(v < 0.0 ? -1 : v > 0.0 ? 1 : 0);
    }
#ifdef USE_MPFR
    if (arg->type == EXPR_MPFR) {
        if (mpfr_zero_p(arg->data.mpfr)) return expr_new_integer(0);
        return expr_new_integer(mpfr_sgn(arg->data.mpfr));
    }
#endif
    int64_t n, d;
    if (is_rational(arg, &n, &d)) {
        /* Canonical Rational has d > 0, but be defensive against an
         * un-canonicalised structural Rational[n, d] reaching here. */
        int s = (n < 0) ^ (d < 0) ? -1 : (n == 0 ? 0 : 1);
        return expr_new_integer(s);
    }

    /* Numeric Complex with real & imag parts: return z / Abs[z]. The
     * evaluator folds the resulting Times[..., Power[Plus[..], -1/2]] to a
     * Complex literal when both parts are concretely numeric. */
    Expr *re, *im;
    if (is_complex(arg, &re, &im) && is_numeric_real(re) && is_numeric_real(im)) {
        /* Zero short-circuit: Sign[0 + 0 I] = 0. */
        bool re_zero = (re->type == EXPR_INTEGER && re->data.integer == 0)
                    || (re->type == EXPR_REAL && re->data.real == 0.0);
        bool im_zero = (im->type == EXPR_INTEGER && im->data.integer == 0)
                    || (im->type == EXPR_REAL && im->data.real == 0.0);
#ifdef USE_MPFR
        if (re->type == EXPR_MPFR) re_zero = mpfr_zero_p(re->data.mpfr) != 0;
        if (im->type == EXPR_MPFR) im_zero = mpfr_zero_p(im->data.mpfr) != 0;
#endif
        if (re_zero && im_zero) return expr_new_integer(0);

#ifdef USE_MPFR
        /* MPFR fast path: avoid the symbolic z * Power[Abs[z], -1] tree
         * (which re-evaluates Abs and Power) and compute the unit-modulus
         * direction directly. */
        if (numeric_any_mpfr(re, im)) {
            long bits = numeric_combined_bits(re, im, 0);
            mpfr_t a_re, a_im, b_re, b_im, abs_v;
            mpfr_init2(a_re, bits); mpfr_init2(a_im, bits);
            mpfr_init2(b_re, bits); mpfr_init2(b_im, bits);
            mpfr_init2(abs_v, bits);
            bool ok_r = get_approx_mpfr(re, a_re, a_im, NULL);
            bool ok_i = get_approx_mpfr(im, b_re, b_im, NULL);
            if (ok_r && ok_i && mpfr_zero_p(a_im) && mpfr_zero_p(b_im)) {
                mpfr_hypot(abs_v, a_re, b_re, MPFR_RNDN);
                /* Above zero short-circuit already handled exact-zero
                 * inputs; abs_v here is strictly positive. */
                Expr* re_out = expr_new_mpfr_bits(bits);
                Expr* im_out = expr_new_mpfr_bits(bits);
                mpfr_div(re_out->data.mpfr, a_re, abs_v, MPFR_RNDN);
                mpfr_div(im_out->data.mpfr, b_re, abs_v, MPFR_RNDN);
                bool im_out_zero = mpfr_zero_p(im_out->data.mpfr);
                mpfr_clear(a_re); mpfr_clear(a_im);
                mpfr_clear(b_re); mpfr_clear(b_im);
                mpfr_clear(abs_v);
                if (im_out_zero) {
                    expr_free(im_out);
                    return re_out;
                }
                return make_complex(re_out, im_out);
            }
            mpfr_clear(a_re); mpfr_clear(a_im);
            mpfr_clear(b_re); mpfr_clear(b_im);
            mpfr_clear(abs_v);
        }
#endif

        Expr* abs_args[1] = { expr_copy(arg) };
        Expr* abs_call = eval_and_free(expr_new_function(expr_new_symbol("Abs"), abs_args, 1));
        Expr* inv_args[2] = { abs_call, expr_new_integer(-1) };
        Expr* inv = eval_and_free(expr_new_function(expr_new_symbol("Power"), inv_args, 2));
        Expr* times_args[2] = { expr_copy(arg), inv };
        return eval_and_free(expr_new_function(expr_new_symbol("Times"), times_args, 2));
    }

    return NULL;
}

Expr* builtin_conjugate(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    if (res->data.function.arg_count != 1) {
        /* Mathematica-compatible argx message; the call is left unevaluated. */
        size_t n = res->data.function.arg_count;
        fprintf(stderr,
                "Conjugate::argx: Conjugate called with %zu argument%s; "
                "1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }
    Expr* arg = res->data.function.args[0];
    /* Conjugate is an involution: Conjugate[Conjugate[z]] -> z. */
    if (head_is(arg, SYM_Conjugate) && arg->data.function.arg_count == 1) {
        Expr* inner = arg->data.function.args[0];
        arg->data.function.args[0] = NULL;
        return inner;
    }
    /* Re, Im, Abs, and Arg are real-valued by construction, so they are
     * Conjugate-fixed regardless of whether their argument numericalizes. */
    if ((head_is(arg, SYM_Re) || head_is(arg, SYM_Im) ||
         head_is(arg, SYM_Abs) || head_is(arg, SYM_Arg)) &&
        arg->data.function.arg_count == 1) {
        return expr_copy(arg);
    }
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        Expr** args_times = malloc(sizeof(Expr*) * 2);
        args_times[0] = expr_new_integer(-1);
        args_times[1] = expr_copy(im);
        Expr* minus_im = expr_new_function(expr_new_symbol("Times"), args_times, 2);
        /* expr_new_function copies the pointers out of args_times into its
         * own internal storage; the caller still owns the temporary array. */
        free(args_times);
        return make_complex(expr_copy(re), minus_im);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_copy(arg);
    }
    /* Symbolic expressions that numericalize to a machine real are
     * Conjugate-fixed (e.g. Sqrt[11], 3/Sqrt[11], Pi, Sqrt[2]*Log[3]).
     * Matches the same reality test used inside complex_decompose. */
    if (is_numeric_real(arg)) {
        return expr_copy(arg);
    }
    Expr *cre, *cim;
    if (complex_decompose(arg, &cre, &cim)) {
        if (is_numeric_real(cre) && is_numeric_real(cim)) {
            Expr* neg_i = make_complex(expr_new_integer(0), expr_new_integer(-1));
            Expr* neg_args[2] = { neg_i, cim };
            Expr* neg_im_term = expr_new_function(expr_new_symbol("Times"), neg_args, 2);
            Expr* plus_args[2] = { cre, neg_im_term };
            Expr* plus = expr_new_function(expr_new_symbol("Plus"), plus_args, 2);
            return eval_and_free(plus);
        }
        expr_free(cre);
        expr_free(cim);
    }
    return NULL;
}

Expr* builtin_arg(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    
    Expr* re = NULL;
    Expr* im = NULL;
    int64_t n, d;
    
#ifdef USE_MPFR
    /* Pure MPFR real: return symbolic 0 or Pi based on sign, preserving
     * exactness rather than coercing to a lossy double. */
    if (arg->type == EXPR_MPFR) {
        int sgn = mpfr_sgn(arg->data.mpfr);
        if (sgn > 0 || mpfr_zero_p(arg->data.mpfr)) return expr_new_integer(0);
        return expr_new_symbol("Pi");
    }
#endif

    if (is_complex(arg, &re, &im)) {
#ifdef USE_MPFR
        /* MPFR-aware Arg: when either component carries MPFR, evaluate
         * via mpfr_atan2 at the working precision rather than the double
         * atan2 fallback below. */
        if (numeric_any_mpfr(re, im)) {
            long bits = numeric_combined_bits(re, im, 0);
            mpfr_t a_re, a_im, b_re, b_im;
            mpfr_init2(a_re, bits); mpfr_init2(a_im, bits);
            mpfr_init2(b_re, bits); mpfr_init2(b_im, bits);
            bool ok_r = get_approx_mpfr(re, a_re, a_im, NULL);
            bool ok_i = get_approx_mpfr(im, b_re, b_im, NULL);
            if (ok_r && ok_i && mpfr_zero_p(a_im) && mpfr_zero_p(b_im)) {
                Expr* result = expr_new_mpfr_bits(bits);
                mpfr_atan2(result->data.mpfr, b_re, a_re, MPFR_RNDN);
                mpfr_clear(a_re); mpfr_clear(a_im);
                mpfr_clear(b_re); mpfr_clear(b_im);
                return result;
            }
            mpfr_clear(a_re); mpfr_clear(a_im);
            mpfr_clear(b_re); mpfr_clear(b_im);
        }
#endif
        // re and im are assigned
    } else if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        re = arg;
    } else {
        return NULL;
    }

    double re_val = 0.0;
    double im_val = 0.0;
    bool inexact = false;

    if (re) {
        if (re->type == EXPR_INTEGER) re_val = (double)re->data.integer;
        else if (re->type == EXPR_REAL) { re_val = re->data.real; inexact = true; }
        else if (is_rational(re, &n, &d)) re_val = (double)n / d;
        else return NULL;
    }
    if (im) {
        if (im->type == EXPR_INTEGER) im_val = (double)im->data.integer;
        else if (im->type == EXPR_REAL) { im_val = im->data.real; inexact = true; }
        else if (is_rational(im, &n, &d)) im_val = (double)n / d;
        else return NULL;
    }
    
    if (re_val == 0.0 && im_val == 0.0) {
        return expr_new_integer(0);
    }
    
    if (!inexact) {
        if (im_val == 0.0) {
            if (re_val > 0.0) return expr_new_integer(0);
            else return expr_new_symbol("Pi");
        }
        if (re_val == 0.0) {
            Expr* args[2];
            args[0] = make_rational(im_val > 0.0 ? 1 : -1, 2);
            args[1] = expr_new_symbol("Pi");
            return expr_new_function(expr_new_symbol("Times"), args, 2);
        }
        if (re_val == im_val) {
            Expr* args[2];
            args[0] = make_rational(re_val > 0.0 ? 1 : -3, 4);
            args[1] = expr_new_symbol("Pi");
            return expr_new_function(expr_new_symbol("Times"), args, 2);
        }
        if (re_val == -im_val) {
            Expr* args[2];
            args[0] = make_rational(re_val > 0.0 ? -1 : 3, 4);
            args[1] = expr_new_symbol("Pi");
            return expr_new_function(expr_new_symbol("Times"), args, 2);
        }

        Expr* args[2];
        args[0] = expr_copy(re);
        args[1] = im ? expr_copy(im) : expr_new_integer(0);
        return expr_new_function(expr_new_symbol("ArcTan"), args, 2);
    }
    
    return expr_new_real(atan2(im_val, re_val));
}
