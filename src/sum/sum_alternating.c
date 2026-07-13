/*
 * sum_alternating.c -- Sum`Alternating: infinite alternating rational sums.
 *
 * Handles  Sum_{k=imin}^Infinity  sigma (-1)^k R(k)  with R a rational function
 * of k and a concrete integer lower bound.  The sign factor may appear in any
 * form (-1)^(a k + b) with a an odd integer; sigma = (-1)^b collects the
 * constant part.  R is partial-fractioned over the rationals and each linear
 * pole term  c0/(b1 k + b0)^m  contributes
 *
 *     c0 b1^{-m} (-1)^imin  LerchPhi[-1, m, imin - rho],   rho = -b0/b1,
 *
 * because  Sum_{k>=imin} (-1)^k/(k + a)^m = (-1)^imin LerchPhi[-1, m, imin + a].
 * The Lerch transcendent Phi(-1, m, .) reduces to elementary constants:
 * Log[2]/Pi at m = 1 (via the digamma / Gauss digamma path) and to Catalan /
 * rational*Pi^m at half-integer arguments (Dirichlet beta).  So, e.g.
 *
 *     Sum[(-1)^(k+1)/k, {k, 1, Infinity}]      -> Log[2]
 *     Sum[(-1)^k/(2k+1), {k, 1, Infinity}]     -> 1/4 (Pi - 4)
 *     Sum[(-1)^k/(2k+1)^2, {k, 0, Infinity}]   -> Catalan
 *
 * Irreducible-quadratic (complex) poles and radical poles are out of scope and
 * left unevaluated (return NULL), never a fabricated value.
 */

#include "sum_internal.h"
#include "symtab.h"
#include "eval.h"
#include "expr.h"
#include "attr.h"
#include "poly.h"
#include "arithmetic.h"   /* is_infinity_sym */
#include "sym_names.h"
#include <stdlib.h>
#include <stdint.h>

/* --- small polynomial helpers (mirror the ones in sum_rational.c) --- */

static bool alt_poly_in(Expr* e, Expr* var) {
    Expr* vars[1] = { var };
    return is_polynomial(e, vars, 1);
}
static int alt_pdeg(Expr* e, Expr* var) {
    Expr* ea[1] = { expr_copy(e) };
    Expr* ex = sum_eval("Expand", ea, 1);
    int d = get_degree_poly(ex, var);
    expr_free(ex);
    return d;
}
/* Power[base, neg int]?  sets *base (alias) and *m = -exp. */
static bool alt_neg_power(Expr* e, Expr** base, int* m) {
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type != EXPR_SYMBOL
        || e->data.function.head->data.symbol != SYM_Power
        || e->data.function.arg_count != 2) return false;
    Expr* ex = e->data.function.args[1];
    if (ex->type != EXPR_INTEGER || ex->data.integer >= 0) return false;
    *base = e->data.function.args[0];
    *m = (int)(-ex->data.integer);
    return true;
}
/* Find the (single) var-dependent pole factor of a partial-fraction term,
 * descending into nested Times.  Returns count; records base/mult in fb, fm. */
static int alt_find_pole(Expr* e, Expr* var, Expr** fb, int* fm) {
    Expr* b; int mm;
    if (alt_neg_power(e, &b, &mm) && !sum_free_of(b, var)) { *fb = b; *fm = mm; return 1; }
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol == SYM_Times) {
        int c = 0;
        for (size_t i = 0; i < e->data.function.arg_count; i++)
            c += alt_find_pole(e->data.function.args[i], var, fb, fm);
        return c;
    }
    return 0;
}

/* c (-1)^imin LerchPhi[-1, m, imin - rho] for one linear pole term c0/base^m,
 * base = b1 var + b0.  Returns owned (unevaluated) Expr*, or NULL if base is not
 * linear in var.  c0 is consumed. */
static Expr* alt_term_contribution(Expr* c0, Expr* base, int m, Expr* var, Expr* imin) {
    if (alt_pdeg(base, var) != 1) { expr_free(c0); return NULL; }

    Expr* zero = sum_int(0);
    Expr* one  = sum_int(1);
    Expr* b0 = sum_subst(base, var, zero);            /* base(0) */
    Expr* b_at1 = sum_subst(base, var, one);          /* base(1) */
    expr_free(zero); expr_free(one);
    Expr* b1 = sum_sub(b_at1, b0);                    /* b1 = base(1) - base(0) */
    expr_free(b_at1);

    /* rho = -b0/b1;  a = imin - rho. */
    Expr* rho = sum_eval("Times",
                    (Expr*[]){ sum_int(-1), b0,
                               expr_new_function(expr_new_symbol(SYM_Power),
                                   (Expr*[]){ expr_copy(b1), sum_int(-1) }, 2) }, 3);
    Expr* a = sum_sub(imin, rho);
    expr_free(rho);

    /* c = c0 b1^{-m}. */
    Expr* c = sum_eval("Times",
                  (Expr*[]){ c0,
                             expr_new_function(expr_new_symbol(SYM_Power),
                                 (Expr*[]){ b1, sum_int(-m) }, 2) }, 2);

    /* (-1)^imin. */
    Expr* sgn = expr_new_function(expr_new_symbol(SYM_Power),
                    (Expr*[]){ sum_int(-1), expr_copy(imin) }, 2);
    /* LerchPhi[-1, m, a]. */
    Expr* lp = expr_new_function(expr_new_symbol(SYM_LerchPhi),
                   (Expr*[]){ sum_int(-1), sum_int(m), a }, 3);

    return expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ c, sgn, lp }, 3);
}

Expr* builtin_sum_alternating(Expr* res);

Expr* builtin_sum_alternating(Expr* res) {
    Expr *f, *var, *imin, *imax;
    bool definite;
    if (!sum_stage_args(res, &f, &var, &imin, &imax, &definite)) return NULL;
    if (!definite) return NULL;
    if (!is_infinity_sym(imax)) return NULL;
    if (imin->type != EXPR_INTEGER) return NULL;

    Expr* fc = expr_copy(f);
    Expr* fn = evaluate(fc);
    expr_free(fc);   /* evaluate() borrows; free the copy we made */

    /* Split fn into the sign exponent(s) (-1)^exp and the rest (R). */
    bool is_times = (fn->type == EXPR_FUNCTION
                     && fn->data.function.head->type == EXPR_SYMBOL
                     && fn->data.function.head->data.symbol == SYM_Times);
    size_t n = is_times ? fn->data.function.arg_count : 1;

    Expr* etot = sum_int(0);      /* sum of the (-1)^. exponents */
    Expr* rprod = sum_int(1);     /* product of the remaining factors = R */
    int signcount = 0;
    for (size_t i = 0; i < n; i++) {
        Expr* g = is_times ? fn->data.function.args[i] : fn;
        Expr* base; int mm;
        bool is_sign = (g->type == EXPR_FUNCTION
            && g->data.function.head->type == EXPR_SYMBOL
            && g->data.function.head->data.symbol == SYM_Power
            && g->data.function.arg_count == 2
            && g->data.function.args[0]->type == EXPR_INTEGER
            && g->data.function.args[0]->data.integer == -1
            && !sum_free_of(g->data.function.args[1], var));
        (void)base; (void)mm;
        if (is_sign) {
            signcount++;
            Expr* t = expr_new_function(expr_new_symbol(SYM_Plus),
                          (Expr*[]){ etot, expr_copy(g->data.function.args[1]) }, 2);
            etot = evaluate(t); expr_free(t);
        } else {
            Expr* t = expr_new_function(expr_new_symbol(SYM_Times),
                          (Expr*[]){ rprod, expr_copy(g) }, 2);
            rprod = evaluate(t); expr_free(t);
        }
    }
    expr_free(fn);
    if (signcount == 0) { expr_free(etot); expr_free(rprod); return NULL; }

    /* etot = a var + b, a an odd integer, b an integer.  sigma = (-1)^b. */
    Expr* acoef = sum_eval("D", (Expr*[]){ expr_copy(etot), expr_copy(var) }, 2);
    Expr* zero = sum_int(0);
    Expr* bconst = sum_subst(etot, var, zero);
    expr_free(zero); expr_free(etot);
    bool ok_sign = (acoef->type == EXPR_INTEGER && (acoef->data.integer % 2 != 0)
                    && bconst->type == EXPR_INTEGER);
    int sigma = 1;
    if (ok_sign) sigma = (bconst->data.integer % 2 == 0) ? 1 : -1;
    expr_free(acoef); expr_free(bconst);
    if (!ok_sign) { expr_free(rprod); return NULL; }

    /* Convergence + rationality gate: R = p/q with deg q >= deg p + 1. */
    Expr* tog = sum_eval("Together", (Expr*[]){ expr_copy(rprod) }, 1);
    Expr* num = sum_eval("Numerator", (Expr*[]){ expr_copy(tog) }, 1);
    Expr* den = sum_eval("Denominator", (Expr*[]){ tog }, 1);
    bool okrat = alt_poly_in(num, var) && alt_poly_in(den, var);
    int degp = okrat ? alt_pdeg(num, var) : -1;
    int degq = okrat ? alt_pdeg(den, var) : -1;
    expr_free(num); expr_free(den);
    if (!okrat || degq < 0 || degp < 0 || degq < degp + 1) {
        expr_free(rprod); return NULL;
    }

    /* Partial-fraction R over Q and sum each linear pole term.
     * (sum_eval consumes rprod.) */
    Expr* ap = sum_eval("Apart", (Expr*[]){ rprod, expr_copy(var) }, 2);
    size_t nt = (ap->type == EXPR_FUNCTION && ap->data.function.head->type == EXPR_SYMBOL
                 && ap->data.function.head->data.symbol == SYM_Plus)
                ? ap->data.function.arg_count : 1;

    Expr** contribs = malloc(sizeof(Expr*) * nt);
    size_t nc = 0;
    bool ok = true;
    for (size_t i = 0; i < nt && ok; i++) {
        Expr* term = (nt > 1) ? ap->data.function.args[i] : ap;
        Expr* base = NULL; int m = 0;
        if (alt_find_pole(term, var, &base, &m) != 1) { ok = false; break; }
        /* residue c0 = Cancel[term * base^m]. */
        Expr* basem = expr_new_function(expr_new_symbol(SYM_Power),
                          (Expr*[]){ expr_copy(base), sum_int(m) }, 2);
        Expr* prod = expr_new_function(expr_new_symbol(SYM_Times),
                         (Expr*[]){ expr_copy(term), basem }, 2);
        Expr* c0 = sum_eval("Cancel", (Expr*[]){ prod }, 1);
        Expr* contrib = alt_term_contribution(c0, base, m, var, imin);
        if (!contrib) { ok = false; break; }
        contribs[nc++] = contrib;
    }
    expr_free(ap);   /* rprod was consumed by Apart above */

    if (!ok) {
        for (size_t i = 0; i < nc; i++) expr_free(contribs[i]);
        free(contribs);
        return NULL;
    }

    Expr* summed = (nc == 1) ? contribs[0]
                             : expr_new_function(expr_new_symbol(SYM_Plus), contribs, nc);
    if (nc != 1) { /* array adopted by Plus */ }
    free(contribs);

    /* result = sigma * summed, then Simplify to collapse the elementary form. */
    Expr* scaled = (sigma == 1) ? summed
        : expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){ sum_int(-1), summed }, 2);
    Expr* ev = evaluate(scaled);
    expr_free(scaled);
    Expr* out = sum_eval("Simplify", (Expr*[]){ ev }, 1);
    return out;
}

void sum_alternating_init(void) {
    symtab_add_builtin("Sum`Alternating", builtin_sum_alternating);
    symtab_get_def("Sum`Alternating")->attributes |= ATTR_PROTECTED | ATTR_READPROTECTED;
    symtab_set_docstring("Sum`Alternating",
        "Sum`Alternating[f, k, imin, Infinity] gives the closed form of an "
        "infinite alternating rational sum sigma (-1)^k R(k) with R rational, via "
        "partial fractions and the Lerch transcendent LerchPhi[-1, m, a]. Results "
        "reduce to elementary constants (Log[2], Pi, Catalan). Returns unevaluated "
        "for non-alternating, non-rational, complex-pole, or divergent inputs.");
}
