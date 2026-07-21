/* test_knowles_tower.c — K0 substrate: the RT_PRIM Liouvillian-primitive tower.
 *
 * Knowles' integrands are transcendental *Liouvillian* — they may already contain
 * li / erf / Ei.  KNOWLES_DESIGN.md §2.1 adds an RT_PRIM tower generator kind so
 * such a kernel becomes a primitive monomial theta with theta' = Dcoef, and closes
 * the tower under each primitive's derivative (Erf -> E^(-u^2), Ei -> E^u,
 * li -> Log[u]) so the lower exp/log kernel is itself a tower monomial.  A merged
 * prim-bearing exponential is split (E^(-x^2 - Erf[x]^2) -> E^(-x^2) E^(-Erf[x]^2))
 * so each factor is a proper monomial.
 *
 * The fundamental correctness property of a differential tower is that the tower
 * derivation COMMUTES with kernel substitution.  We verify it in the ORIGINAL
 * kernel domain (back-substituting the tower variables to their kernels, so a power
 * like t0^2 becomes E^(-2x^2) and cancels naturally):
 *
 *     backsubst( D_tower[ subst(f) ] )  ==  D[f, x] .
 *
 * This exercises every generator's Dcoef, including RT_PRIM.  We also assert an
 * RT_PRIM generator was actually created for each SF kernel, and that purely
 * elementary towers are unchanged (0 prim, byte-identical to before).
 */

#include "eval.h"
#include "parse.h"
#include "expr.h"
#include "symtab.h"
#include "core.h"
#include "print.h"
#include "risch_tower.h"
#include "test_utils.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdbool.h>

static int failures = 0;

static size_t count_prim(const RtTower* T) {
    size_t k = 0;
    for (size_t i = 0; i < T->n; i++) if (T->kind[i] == RT_PRIM) k++;
    return k;
}

/* Build the tower of f (a WL string) over x; verify the commutation identity and
 * that it has exactly `want_prim` RT_PRIM generators and depth >= want_n. */
static void check_tower(const char* fstr, size_t want_prim, size_t want_n) {
    Expr* x = expr_new_symbol("x");
    Expr* f0 = parse_expression(fstr);
    Expr* fe = evaluate(f0); expr_free(f0);
    /* Pre-normalise like the integrator: split prim-bearing exponentials so the
     * integrand's kernels match the tower's. */
    Expr* fx = rt_expand_exp_sums(fe);
    expr_free(fe);

    RtTower T;
    bool built = rt_tower_build_min(fx, x, &T, 1);
    if (!built) {
        printf("  FAIL [%s]: tower did not build\n", fstr);
        failures++;
        rt_tower_free(&T); expr_free(fx); expr_free(x);
        return;
    }

    size_t np = count_prim(&T);
    if (np != want_prim) {
        printf("  FAIL [%s]: RT_PRIM count %zu, expected %zu (depth %zu)\n",
               fstr, np, want_prim, T.n);
        failures++;
    }
    if (T.n < want_n) {
        printf("  FAIL [%s]: tower depth %zu < expected %zu\n", fstr, T.n, want_n);
        failures++;
    }

    /* lhs = D_tower[subst(fx)] in tower variables. */
    Expr* fsub = rt_subst_kernels(fx, &T);
    Expr* dtow = rt_tower_deriv(fsub, &T, x);
    expr_free(fsub);

    /* Back-substitute t_i -> kernel_i simultaneously (kernels carry no tower var). */
    Expr** rules = malloc((T.n ? T.n : 1) * sizeof(Expr*));
    for (size_t i = 0; i < T.n; i++)
        rules[i] = expr_new_function(expr_new_symbol("Rule"),
            (Expr*[]){ expr_copy(T.t[i]), expr_copy(T.kernel[i]) }, 2);
    Expr* rlist = expr_new_function(expr_new_symbol("List"), rules, T.n);
    free(rules);
    Expr* back = evaluate(expr_new_function(expr_new_symbol("ReplaceAll"),
        (Expr*[]){ dtow, rlist }, 2));           /* adopts dtow, rlist */

    Expr* df = evaluate(expr_new_function(expr_new_symbol("D"),
        (Expr*[]){ expr_copy(fx), expr_copy(x) }, 2));
    Expr* diff = expr_new_function(expr_new_symbol("Subtract"),
        (Expr*[]){ back, df }, 2);
    Expr* sm = evaluate(expr_new_function(expr_new_symbol("Simplify"),
        (Expr*[]){ diff }, 1));                   /* adopts diff */
    char* s = expr_to_string_fullform(sm);
    bool zero = (strcmp(s, "0") == 0);
    if (!zero) { printf("  FAIL [%s]: commutation nonzero -> %s\n", fstr, s); failures++; }
    else printf("  ok   [%s]  (depth %zu, %zu prim)\n", fstr, T.n, np);
    free(s); expr_free(sm);

    rt_tower_free(&T); expr_free(fx); expr_free(x);
}

int main(void) {
    core_init();
    printf("Running test: knowles_tower (K0 RT_PRIM substrate)\n");

    /* --- Primitive SF kernels build and derive correctly --------------------- */
    check_tower("Erf[x]", 1, 2);                 /* E^(-x^2) < Erf[x]          */
    check_tower("Erfi[x]", 1, 2);                /* E^(x^2)  < Erfi[x]         */
    check_tower("Erfc[x]", 1, 2);                /* E^(-x^2) < Erfc[x]         */
    check_tower("ExpIntegralEi[1/x]", 1, 2);     /* E^(1/x)  < Ei[1/x]         */
    check_tower("LogIntegral[x]", 1, 2);         /* Log[x]   < li[x]           */

    /* --- Nested & mixed Liouvillian towers ---------------------------------- */
    /* Erf[Erf[x]]: E^(-x^2) < Erf[x] < E^(-Erf[x]^2) < Erf[Erf[x]]  (2 prim). */
    check_tower("Erf[Erf[x]]", 2, 4);
    /* li[li[x]] (1986 example answer kernel): Log[x] < li[x] < Log[li[x]] < li[li[x]]. */
    check_tower("LogIntegral[LogIntegral[x]]", 2, 4);
    /* Part II Ex 4.3 integrand shape: Exp[-x^2] Erf[x]. */
    check_tower("Exp[-x^2] Erf[x]", 1, 2);
    /* Part II Ex 4.1 integrand kernel: Exp[-x^2 - Erf[x]^2] (split into 3 monomials). */
    check_tower("Exp[-x^2 - Erf[x]^2]", 1, 3);

    /* --- Non-regression: purely elementary towers are unchanged (0 prim) ----- */
    check_tower("Log[x]", 0, 1);
    check_tower("Exp[x]", 0, 1);
    check_tower("Exp[x] Log[x]", 0, 2);
    check_tower("Log[Log[x]]", 0, 2);

    if (failures) { printf("knowles_tower: %d FAIL\n", failures); return 1; }
    printf("knowles_tower: all passed\n");
    return 0;
}
