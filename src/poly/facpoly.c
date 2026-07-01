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

void facpoly_init(void) {
    symtab_add_builtin("FactorSquareFree", builtin_factorsquarefree);
    symtab_get_def("FactorSquareFree")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("Factor", builtin_factor);
    symtab_get_def("Factor")->attributes |= ATTR_PROTECTED | ATTR_LISTABLE;
    symtab_add_builtin("FactorTerms", builtin_factorterms);
    symtab_get_def("FactorTerms")->attributes |= ATTR_PROTECTED;
    symtab_add_builtin("FactorTermsList", builtin_factortermslist);
    symtab_get_def("FactorTermsList")->attributes |= ATTR_PROTECTED;
}
