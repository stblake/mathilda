/*
 * sum_gosper.c -- Sum`Gosper: Gosper's indefinite hypergeometric summation,
 * plus the DifferenceDelta forward-difference operator.
 *
 * Given a hypergeometric term t(i) (one whose term ratio t(i+1)/t(i) is a
 * rational function of i), Gosper's algorithm finds a hypergeometric
 * antidifference F with F(i+1)-F(i) = t(i), or proves none exists.  The output
 * has the shape F = R(i) t(i) with R rational, so no new special functions are
 * needed.
 *
 *   1. r(i) = t(i+1)/t(i); require it rational (Simplify reduces factorial
 *      ratios, then Together gives num/den polynomials a, b).
 *   2. Gosper-Petkovsek normal form r = (a/b)(c(i+1)/c(i)) with
 *      gcd(a(i), b(i+h)) = 1 for all integers h >= 0, via the dispersion set
 *      (h with deg gcd(a(i), b(i+h)) > 0) and gcd peeling.
 *   3. Solve a(i) x(i+1) - b(i-1) x(i) = c(i) for a polynomial x by undetermined
 *      coefficients (SolveAlways).  No solution => t is not Gosper-summable.
 *   4. Antidifference F(i) = (b(i-1)/c(i)) x(i) t(i).
 *
 *   Sum`Gosper[f, i]              -> F(i)                 (indefinite)
 *   Sum`Gosper[f, i, imin, imax]  -> F(imax+1) - F(imin)  (definite, finite)
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "sym_names.h"
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>

#define GOSPER_DISPERSION_MAX 64   /* cap on dispersion-set search */

static int gosper_counter = 0;

/* ---- small polynomial helpers built on existing builtins ---- */

/* degree of e in var (expanding first). */
static int pdeg(Expr* e, Expr* var) {
    Expr* ea[1] = { expr_copy(e) };
    Expr* ex = sum_eval("Expand", ea, 1);
    int d = get_degree_poly(ex, var);
    expr_free(ex);
    return d;
}

/* e is a polynomial in var? */
static bool poly_in(Expr* e, Expr* var) {
    Expr* vars[1] = { var };
    return is_polynomial(e, vars, 1);
}

/* PolynomialGCD[a, b] */
static Expr* pgcd(Expr* a, Expr* b) {
    Expr* args[2] = { expr_copy(a), expr_copy(b) };
    return sum_eval("PolynomialGCD", args, 2);
}

/* PolynomialQuotient[a, b, var] */
static Expr* pquot(Expr* a, Expr* b, Expr* var) {
    Expr* args[3] = { expr_copy(a), expr_copy(b), expr_copy(var) };
    return sum_eval("PolynomialQuotient", args, 3);
}

/* e /. var -> var + k */
static Expr* shift_var(Expr* e, Expr* var, int k) {
    Expr* nv = expr_new_function(expr_new_symbol("Plus"),
                   (Expr*[]){ expr_copy(var), sum_int(k) }, 2);
    Expr* r = sum_subst(e, var, nv);
    expr_free(nv);
    return r;
}

/*
 * Gosper's algorithm.  Returns the indefinite antidifference F(var) such that
 * F(var+1) - F(var) = t(var), or NULL if t is not a hypergeometric term or not
 * Gosper-summable.
 */
static Expr* gosper_antidiff(Expr* t, Expr* var) {
    /* 1. term ratio r = t(var+1)/t(var), reduced to a rational function. */
    Expr* tshift = shift_var(t, var, 1);
    Expr* ratio_raw = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ tshift,
                                     expr_new_function(expr_new_symbol("Power"),
                                         (Expr*[]){ expr_copy(t), sum_int(-1) }, 2) }, 2);
    Expr* sarg[1] = { ratio_raw };
    Expr* simp = sum_eval("Simplify", sarg, 1);
    Expr* targ[1] = { simp };
    Expr* ratio = sum_eval("Together", targ, 1);

    Expr* na[1] = { expr_copy(ratio) };
    Expr* num = sum_eval("Numerator", na, 1);
    Expr* dna[1] = { expr_copy(ratio) };
    Expr* den = sum_eval("Denominator", dna, 1);
    expr_free(ratio);

    if (!poly_in(num, var) || !poly_in(den, var)) {
        expr_free(num); expr_free(den);
        return NULL;
    }

    /* 2. Gosper-Petkovsek normal form.  a, b, c polynomials. */
    Expr* a = num;            /* takes ownership */
    Expr* b = den;
    Expr* c = sum_int(1);
    Expr* f0 = expr_copy(a);
    Expr* g0 = expr_copy(b);

    int Z[GOSPER_DISPERSION_MAX];
    int zc = 0;
    for (int h = 1; h <= GOSPER_DISPERSION_MAX; h++) {
        Expr* g0h = shift_var(g0, var, h);
        Expr* s = pgcd(f0, g0h);
        if (pdeg(s, var) > 0) Z[zc++] = h;
        expr_free(g0h); expr_free(s);
    }
    expr_free(f0); expr_free(g0);

    for (int zi = 0; zi < zc; zi++) {
        int h = Z[zi];
        Expr* bh = shift_var(b, var, h);
        Expr* s = pgcd(a, bh);
        expr_free(bh);
        if (pdeg(s, var) <= 0) { expr_free(s); continue; }

        Expr* a2 = pquot(a, s, var);
        expr_free(a); a = a2;

        Expr* sh = shift_var(s, var, -h);
        Expr* b2 = pquot(b, sh, var);
        expr_free(sh);
        expr_free(b); b = b2;

        for (int i = 1; i <= h; i++) {
            Expr* si = shift_var(s, var, -i);
            Expr* tc = expr_new_function(expr_new_symbol("Times"),
                           (Expr*[]){ c, si }, 2);
            c = evaluate(tc);
            expr_free(tc);
        }
        expr_free(s);
    }

    /* 3. Solve a(var) x(var+1) - b(var-1) x(var) = c(var) for polynomial x. */
    Expr* bb = shift_var(b, var, -1);          /* b(var-1) */
    int da = pdeg(a, var), db = pdeg(b, var), dc = pdeg(c, var);
    int Dmax = dc + (da > db ? da : db) + 2;
    if (Dmax < 0) Dmax = 0;

    Expr* x_sol = NULL;
    int id = ++gosper_counter;

    for (int d = 0; d <= Dmax && !x_sol; d++) {
        /* x = sum_{k=0}^{d} u_k var^k with fresh coefficient symbols. */
        char name[64];
        Expr** csym = malloc(sizeof(Expr*) * (d + 1));
        Expr** xt = malloc(sizeof(Expr*) * (d + 1));
        for (int k = 0; k <= d; k++) {
            snprintf(name, sizeof name, "Sum`gx`%d`%d", id, k);
            csym[k] = expr_new_symbol(name);
            Expr* pw = expr_new_function(expr_new_symbol("Power"),
                           (Expr*[]){ expr_copy(var), sum_int(k) }, 2);
            xt[k] = expr_new_function(expr_new_symbol("Times"),
                        (Expr*[]){ expr_copy(csym[k]), pw }, 2);
        }
        Expr** xa = malloc(sizeof(Expr*) * (d + 1));
        for (int k = 0; k <= d; k++) xa[k] = expr_copy(xt[k]);
        Expr* x = expr_new_function(expr_new_symbol("Plus"), xa, d + 1);
        free(xa);
        for (int k = 0; k <= d; k++) expr_free(xt[k]);
        free(xt);

        /* eqn = a*x(var+1) - bb*x - c */
        Expr* xshift = shift_var(x, var, 1);
        Expr* term1 = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ expr_copy(a), xshift }, 2);
        Expr* term2 = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ sum_int(-1), expr_copy(bb), expr_copy(x) }, 3);
        Expr* term3 = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ sum_int(-1), expr_copy(c) }, 2);
        Expr* eqraw = expr_new_function(expr_new_symbol("Plus"),
                          (Expr*[]){ term1, term2, term3 }, 3);
        Expr* eea[1] = { eqraw };
        Expr* eqn = sum_eval("Expand", eea, 1);

        Expr* eq = expr_new_function(expr_new_symbol("Equal"),
                       (Expr*[]){ eqn, sum_int(0) }, 2);
        Expr* saa[2] = { eq, expr_copy(var) };
        Expr* sol = sum_eval("SolveAlways", saa, 2);

        bool solved = false;
        if (sol->type == EXPR_FUNCTION
            && sol->data.function.head->type == EXPR_SYMBOL
            && sol->data.function.head->data.symbol == SYM_List
            && sol->data.function.arg_count >= 1) {
            Expr* rules = sol->data.function.args[0];
            if (rules->type == EXPR_FUNCTION
                && rules->data.function.arg_count >= 1) {
                /* apply the solution, then zero any free coefficients. */
                Expr* ra[2] = { expr_copy(x), expr_copy(rules) };
                Expr* xr = sum_eval("ReplaceAll", ra, 2);
                /* zero leftovers */
                Expr** zr = malloc(sizeof(Expr*) * (d + 1));
                for (int k = 0; k <= d; k++)
                    zr[k] = expr_new_function(expr_new_symbol("Rule"),
                                (Expr*[]){ expr_copy(csym[k]), sum_int(0) }, 2);
                Expr* zrl = expr_new_function(expr_new_symbol("List"), zr, d + 1);
                free(zr);
                Expr* za[2] = { xr, zrl };
                x_sol = sum_eval("ReplaceAll", za, 2);
                solved = true;
            }
        }
        expr_free(sol);
        expr_free(x);
        for (int k = 0; k <= d; k++) expr_free(csym[k]);
        free(csym);
        (void)solved;
    }

    Expr* result = NULL;
    if (x_sol) {
        /* 4. F = (b(var-1)/c) * x * t, cleaned with Cancel. */
        Expr* Rcoef = expr_new_function(expr_new_symbol("Times"),
                          (Expr*[]){ expr_copy(bb), x_sol,
                                     expr_new_function(expr_new_symbol("Power"),
                                         (Expr*[]){ expr_copy(c), sum_int(-1) }, 2) }, 3);
        Expr* Fraw = expr_new_function(expr_new_symbol("Times"),
                         (Expr*[]){ Rcoef, expr_copy(t) }, 2);
        Expr* F = evaluate(Fraw);
        expr_free(Fraw);
        Expr* ca[1] = { F };
        result = sum_eval("Cancel", ca, 1);
    }

    expr_free(bb);
    expr_free(a); expr_free(b); expr_free(c);
    return result;
}

Expr* builtin_sum_gosper(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;

    if (definite && imax->type == EXPR_SYMBOL && imax->data.symbol == SYM_Infinity)
        return NULL;

    Expr* F = gosper_antidiff(f, var);
    if (!F) return NULL;

    if (!definite) return F;

    /* Definite: F(imax+1) - F(imin), cleaned with Cancel. */
    Expr* up = expr_new_function(expr_new_symbol("Plus"),
                   (Expr*[]){ expr_copy(imax), sum_int(1) }, 2);
    Expr* Fhi = sum_subst(F, var, up);
    Expr* Flo = sum_subst(F, var, imin);
    expr_free(up);
    expr_free(F);
    Expr* diff = sum_sub(Fhi, Flo);
    expr_free(Fhi); expr_free(Flo);
    Expr* ca[1] = { diff };
    return sum_eval("Cancel", ca, 1);
}

/* ---- DifferenceDelta[f, i] = (f /. i -> i+1) - f ---- */
Expr* builtin_differencedelta(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count != 2)
        return NULL;
    Expr* f   = res->data.function.args[0];
    Expr* var = res->data.function.args[1];
    if (var->type != EXPR_SYMBOL) return NULL;

    Expr* shifted = shift_var(f, var, 1);
    Expr* diff = sum_sub(shifted, f);
    expr_free(shifted);
    Expr* ea[1] = { diff };
    return sum_eval("Expand", ea, 1);
}

void sum_gosper_init(void) {
    symtab_add_builtin("Sum`Gosper", builtin_sum_gosper);
    symtab_get_def("Sum`Gosper")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`Gosper",
        "Sum`Gosper[f, i] gives the indefinite sum of a hypergeometric term f "
        "in i via Gosper's algorithm; Sum`Gosper[f, i, imin, imax] gives the "
        "finite definite sum. Returns unevaluated if f is not a hypergeometric "
        "term or is not Gosper-summable.");

    symtab_add_builtin("DifferenceDelta", builtin_differencedelta);
    symtab_get_def("DifferenceDelta")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("DifferenceDelta",
        "DifferenceDelta[f, i] gives the forward difference (f /. i -> i+1) - f, "
        "the discrete analogue of D. It is the left inverse of indefinite Sum.");
}
