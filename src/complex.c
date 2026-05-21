#include "complex.h"
#include "symtab.h"
#include "eval.h"
#include "arithmetic.h"
#include "numeric.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

static bool head_is(Expr* e, const char* name) {
    return e && e->type == EXPR_FUNCTION &&
           e->data.function.head &&
           e->data.function.head->type == EXPR_SYMBOL &&
           strcmp(e->data.function.head->data.symbol, name) == 0;
}

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

    if (head_is(e, "Plus")) {
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

    if (head_is(e, "Times")) {
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
    symtab_add_builtin("Conjugate", builtin_conjugate);
    symtab_add_builtin("Arg", builtin_arg);

    symtab_get_def("Re")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Im")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("ReIm")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Abs")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Conjugate")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
    symtab_get_def("Arg")->attributes |= (ATTR_PROTECTED | ATTR_NUMERICFUNCTION | ATTR_LISTABLE);
}

Expr* builtin_re(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        return expr_copy(re);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_copy(arg);
    }
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
    Expr *re, *im;
    if (is_complex(arg, &re, &im)) {
        return expr_copy(im);
    }
    int64_t n, d;
    if (arg->type == EXPR_INTEGER || arg->type == EXPR_REAL || is_rational(arg, &n, &d)) {
        return expr_new_integer(0);
    }
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
    if (is_rational(arg, &n, &d)) {
        return make_rational(n < 0 ? -n : n, d);
    }
    return NULL;
}

Expr* builtin_conjugate(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 1) return NULL;
    Expr* arg = res->data.function.args[0];
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
    
    if (is_complex(arg, &re, &im)) {
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
