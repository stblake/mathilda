#include "facpoly.h"
#include "poly.h"
#include "eval.h"
#include "expand.h"
#include "symtab.h"
#include "attr.h"
#include "arithmetic.h"
#include "internal.h"
#include "parse.h"
#include "rationalize.h"
#include "zupoly.h"
#include "bpoly.h"
#include "mpoly.h"
#include "mvfactor.h"
#include "mvfactor3.h"
#include "sym_names.h"
#include "qafactor.h"
#include "flint_bridge.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <limits.h>
#include <gmp.h>

/* facpoly.c is split into a set of .inc fragments included into this single
 * translation unit.  All in-fragment statics keep file scope; the forward
 * declarations below let the fragments be #included in any order.
 *
 *   facpoly_squarefree.inc        FactorSquareFree pipeline.
 *   facpoly_simple.inc            int_root, binomial, degree-one heuristics.
 *   facpoly_monomial_content.inc  Variable-monomial GCD extraction.
 *   facpoly_irreducibility.inc    Image-based irreducibility probe + factor_roots.
 *   facpoly_bv_hensel.inc         Bivariate Hensel lift (Stages 1-3).
 *   facpoly_z_split.inc           n-variate specialise-and-divide.
 *   facpoly_tv_hensel.inc         Trivariate Hensel via MPoly.
 *   facpoly_heuristic.inc         heuristic_factor recursive dispatcher.
 *   facpoly_memo.inc              Per-Simplify Factor result cache.
 *   facpoly_factor_builtin.inc    builtin_factor entry point.
 *   facpoly_factorterms.inc       FactorTerms / FactorTermsList.
 *   facpoly_bz_uni.inc            Univariate Berlekamp-Zassenhaus pipeline.
 */
static Expr* heuristic_factor(Expr* P);
static Expr* factor_monomial_content(Expr* P);
static Expr* factor_binomial(Expr* P);
static Expr* factor_degree_one(Expr* P, Expr** vars, size_t v_count);
static Expr* factor_bivariate_via_hensel(Expr* P, Expr** vars);
static Expr* factor_via_z_independent_split(Expr* P, Expr** vars, size_t v_count);
static Expr* factor_trivariate_via_mhensel(Expr* P, Expr** vars, size_t v_count);
static Expr* factor_roots(Expr* P, Expr** vars, size_t v_count);
static bool is_likely_irreducible_multivariate(Expr* P, Expr** vars, size_t v_count);
static Expr* eval_others_at_alpha(Expr* P, Expr** vars, size_t v_count,
                                  size_t main_idx, int64_t alpha);
static bool univariate_squarefree(Expr* u, Expr* var);
static FactorMemo* factor_memo_top(void);
static const Expr* factor_memo_get(FactorMemo* m, Expr* key);
static void factor_memo_put(FactorMemo* m, Expr* key, Expr* value);

Expr* bz_factor_to_expr(Expr* P, Expr* var);

#include "facpoly_squarefree.inc"
#include "facpoly_simple.inc"
#include "facpoly_monomial_content.inc"
#include "facpoly_irreducibility.inc"
#include "facpoly_bv_hensel.inc"
#include "facpoly_z_split.inc"
#include "facpoly_tv_hensel.inc"
#include "facpoly_heuristic.inc"
#include "facpoly_memo.inc"
#include "facpoly_factor_builtin.inc"
#include "facpoly_factorterms.inc"
#include "facpoly_bz_uni.inc"

/* ===================================================================== */
/* PolynomialSqrt[p] / PolynomialSqrt[p, x] — the polynomial square root.  */
/* ===================================================================== */
/* Returns s with s^2 == p exactly when p is a perfect square (every
 * non-constant irreducible factor occurs to an even multiplicity), else the
 * sentinel `$Failed`.  The numeric content is carried through Sqrt (so
 * PolynomialSqrt[2(1+x)^2] = Sqrt[2] (1+x) over the reals / F-bar), which keeps
 * the routine total over an algebraically closed constant field — exactly the
 * setting Cherry's completing-square Erf-argument step needs (r_i = Sqrt[p+beta q]).
 * Correct by construction: the candidate s is accepted only after the exact
 * Expand[s^2 - p] == 0 check, so a mis-factorisation can never return a wrong root. */

static bool ps_is_numeric(const Expr* e) {
    if (!e) return false;
    switch (e->type) {
        case EXPR_INTEGER: case EXPR_BIGINT: case EXPR_REAL: return true;
        case EXPR_FUNCTION:
            return e->data.function.head->type == EXPR_SYMBOL &&
                   (e->data.function.head->data.symbol.name == SYM_Rational ||
                    e->data.function.head->data.symbol.name == SYM_Complex);
        default: return false;
    }
}

Expr* builtin_polynomialsqrt(Expr* res) {
    if (res->type != EXPR_FUNCTION || res->data.function.arg_count < 1
        || res->data.function.arg_count > 2)
        return NULL;
    Expr* p = res->data.function.args[0];

    /* 0 -> 0. */
    Expr* pe = expr_expand(p);
    if (is_zero_poly(pe)) { expr_free(pe); return expr_new_integer(0); }
    expr_free(pe);

    /* Factor into  constant * prod(base^exp). */
    Expr** fa = malloc(sizeof(Expr*)); fa[0] = expr_copy(p);
    Expr* f = internal_factor(fa, 1); free(fa);
    if (!f) return NULL;

    Expr** items; size_t nit; Expr* single[1];
    if (f->type == EXPR_FUNCTION && f->data.function.head->type == EXPR_SYMBOL
        && f->data.function.head->data.symbol.name == SYM_Times) {
        items = f->data.function.args; nit = f->data.function.arg_count;
    } else { single[0] = f; items = single; nit = 1; }

    Expr* konst = expr_new_integer(1);                       /* numeric content */
    Expr** halves = malloc((nit ? nit : 1) * sizeof(Expr*)); size_t nh = 0;
    bool ok = true;
    for (size_t i = 0; i < nit && ok; i++) {
        Expr* t = items[i];
        if (ps_is_numeric(t)) {
            konst = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                (Expr*[]){ konst, expr_copy(t) }, 2));
        } else if (t->type == EXPR_FUNCTION
                   && t->data.function.head->type == EXPR_SYMBOL
                   && t->data.function.head->data.symbol.name == SYM_Power
                   && t->data.function.arg_count == 2
                   && t->data.function.args[1]->type == EXPR_INTEGER) {
            int64_t e = t->data.function.args[1]->data.integer;
            Expr* base = t->data.function.args[0];
            if (ps_is_numeric(base)) {
                konst = eval_and_free(expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ konst, expr_copy(t) }, 2));
            } else if (e % 2 == 0) {
                halves[nh++] = expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(base), expr_new_integer(e / 2) }, 2);
            } else {
                ok = false;                       /* odd multiplicity, non-constant */
            }
        } else {
            ok = false;                           /* bare irreducible (multiplicity 1) */
        }
    }
    expr_free(f);

    Expr* result = NULL;
    if (ok) {
        Expr* sk = eval_and_free(expr_new_function(expr_new_symbol("Power"),
            (Expr*[]){ konst, expr_new_function(expr_new_symbol("Rational"),
                (Expr*[]){ expr_new_integer(1), expr_new_integer(2) }, 2) }, 2));
        Expr** facs = malloc((nh + 1) * sizeof(Expr*)); size_t nf = 0;
        facs[nf++] = sk;
        for (size_t i = 0; i < nh; i++) facs[nf++] = halves[i];
        Expr* s = (nf == 1) ? facs[0]
                : eval_and_free(expr_new_function(expr_new_symbol("Times"), facs, nf));
        free(facs);
        /* Exact certificate: Expand[s^2 - p] == 0. */
        Expr* chk = eval_and_free(expr_new_function(expr_new_symbol("Expand"),
            (Expr*[]){ expr_new_function(expr_new_symbol("Plus"),
                (Expr*[]){ expr_new_function(expr_new_symbol("Power"),
                    (Expr*[]){ expr_copy(s), expr_new_integer(2) }, 2),
                  expr_new_function(expr_new_symbol("Times"),
                    (Expr*[]){ expr_new_integer(-1), expr_copy(p) }, 2) }, 2) }, 1));
        if (chk && is_zero_poly(chk)) result = s; else expr_free(s);
        if (chk) expr_free(chk);
    } else {
        expr_free(konst);
        for (size_t i = 0; i < nh; i++) expr_free(halves[i]);
    }
    free(halves);
    return result ? result : expr_new_symbol("$Failed");
}

void facpoly_init(void) {
    symtab_add_builtin("FactorSquareFree", builtin_factorsquarefree);
    symtab_get_def("FactorSquareFree")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("Factor", builtin_factor);
    symtab_get_def("Factor")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("FactorTerms", builtin_factorterms);
    symtab_get_def("FactorTerms")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FactorTermsList", builtin_factortermslist);
    symtab_get_def("FactorTermsList")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("PolynomialSqrt", builtin_polynomialsqrt);
    symtab_get_def("PolynomialSqrt")->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("PolynomialSqrt",
        "PolynomialSqrt[p] gives a polynomial s with s^2 == p when p is a perfect "
        "square (every non-constant irreducible factor has even multiplicity; the "
        "numeric content is carried through Sqrt), and $Failed otherwise. "
        "PolynomialSqrt[p, x] treats p as a polynomial in x.");
}
