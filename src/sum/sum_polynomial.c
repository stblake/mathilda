/*
 * sum_polynomial.c -- Sum`Polynomial: closed-form summation of polynomials.
 *
 * For a polynomial f(i) of degree d we build the antidifference F (the unique
 * polynomial with F(i+1)-F(i) = f(i) and F(0) = 0) using Newton's forward
 * difference formula in the falling-factorial basis -- no Bernoulli numbers
 * required:
 *
 *     f(i) = sum_{k=0}^{d}  Delta^k f(0) * Binomial[i, k]
 *     F(i) = sum_{k=0}^{d}  Delta^k f(0) * Binomial[i, k+1]
 *          = sum_{k=0}^{d} (Delta^k f(0) / (k+1)!) * i(i-1)...(i-k)
 *
 * because the discrete antidifference of the falling factorial i^{(k)} is
 * i^{(k+1)}/(k+1).  The forward differences Delta^k f(0) are computed from the
 * d+1 sample values f(0..d) by repeated differencing, so symbolic coefficients
 * flow through unchanged.
 *
 *   Sum`Polynomial[f, i]              -> F(i)                 (indefinite)
 *   Sum`Polynomial[f, i, imin, imax]  -> F(imax+1) - F(imin)  (definite)
 *
 * The result is run through Factor so that, e.g., Sum[i^2,i] prints as the
 * familiar 1/6 (-1+i) i (-1+2 i).
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "sym_names.h"
#include <stdlib.h>
#include <stdint.h>

/* Build the falling-factorial product i(i-1)...(i-k) = prod_{j=0}^{k}(i-j). */
static Expr* falling_factorial(Expr* var, int k) {
    Expr** factors = malloc(sizeof(Expr*) * (k + 1));
    for (int j = 0; j <= k; j++) {
        if (j == 0) {
            factors[j] = expr_copy(var);
        } else {
            factors[j] = expr_new_function(expr_new_symbol(SYM_Plus),
                (Expr*[]){ expr_copy(var), sum_int(-j) }, 2);
        }
    }
    Expr* prod = expr_new_function(expr_new_symbol(SYM_Times), factors, k + 1);
    free(factors);
    return prod;
}

/* Build the antidifference polynomial F(var) with F(0)=0, Delta F = f.
 * Returns NULL if f is not a polynomial in var. */
static Expr* antidifference(Expr* f_in, Expr* var) {
    Expr* vars[1] = { var };
    if (!is_polynomial(f_in, vars, 1)) return NULL;

    /* Expand so the degree is read correctly for forms like (i+3)^5. */
    Expr* fe[1] = { expr_copy(f_in) };
    Expr* f = sum_eval("Expand", fe, 1);

    int d = get_degree_poly(f, var);
    if (d < 0) d = 0;

    /* Sample f at var = 0, 1, ..., d (symbolic coefficients survive). */
    Expr** vals = malloc(sizeof(Expr*) * (d + 1));
    for (int m = 0; m <= d; m++) {
        Expr* mm = sum_int(m);
        vals[m] = sum_subst(f, var, mm);
        expr_free(mm);
    }
    expr_free(f);

    /* Forward-difference table in place: after step k, vals[0] = Delta^k f(0).
     * delta0[k] captures the leading difference at each level. */
    Expr** delta0 = malloc(sizeof(Expr*) * (d + 1));
    for (int k = 0; k <= d; k++) {
        delta0[k] = expr_copy(vals[0]);
        /* vals[m] <- vals[m+1] - vals[m] for m = 0 .. d-k-1 */
        for (int m = 0; m + k < d; m++) {
            Expr* diff = sum_sub(vals[m + 1], vals[m]);
            expr_free(vals[m]);
            vals[m] = diff;
        }
    }
    for (int m = 0; m <= d; m++) expr_free(vals[m]);
    free(vals);

    /* F = sum_k delta0[k] * fallingfactorial(k) / (k+1)! */
    Expr** terms = malloc(sizeof(Expr*) * (d + 1));
    for (int k = 0; k <= d; k++) {
        int64_t fact = 1;
        for (int j = 2; j <= k + 1; j++) fact *= j;   /* (k+1)! */
        Expr* prod = falling_factorial(var, k);
        Expr* invf = expr_new_function(expr_new_symbol(SYM_Power),
                        (Expr*[]){ sum_int(fact), sum_int(-1) }, 2);
        terms[k] = expr_new_function(expr_new_symbol(SYM_Times),
                        (Expr*[]){ delta0[k], prod, invf }, 3);
    }
    free(delta0);

    Expr* sum = expr_new_function(expr_new_symbol(SYM_Plus), terms, d + 1);
    free(terms);
    Expr* F = evaluate(sum);
    expr_free(sum);
    return F;
}

Expr* builtin_sum_polynomial(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;

    /* A polynomial summed to Infinity diverges; fall through (stays held). */
    if (definite && imax->type == EXPR_SYMBOL && imax->data.symbol.name == SYM_Infinity)
        return NULL;

    Expr* F = antidifference(f, var);
    if (!F) return NULL;

    if (!definite) {
        Expr* r = sum_factor(F);
        expr_free(F);
        return r;
    }

    /* Definite: F(imax+1) - F(imin). */
    Expr* up = expr_new_function(expr_new_symbol(SYM_Plus),
                  (Expr*[]){ expr_copy(imax), sum_int(1) }, 2);
    Expr* Fhi = sum_subst(F, var, up);
    Expr* Flo = sum_subst(F, var, imin);
    expr_free(up);
    expr_free(F);

    Expr* diff = sum_sub(Fhi, Flo);
    expr_free(Fhi);
    expr_free(Flo);

    Expr* r = sum_factor(diff);
    expr_free(diff);
    return r;
}

void sum_polynomial_init(void) {
    symtab_add_builtin("Sum`Polynomial", builtin_sum_polynomial);
    symtab_get_def("Sum`Polynomial")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`Polynomial",
        "Sum`Polynomial[f, i] gives the indefinite sum of a polynomial f in i "
        "(its antidifference). Sum`Polynomial[f, i, imin, imax] gives the "
        "definite sum. Uses Newton forward differences in the falling-factorial "
        "basis; returns unevaluated if f is not a polynomial in i.");
}
