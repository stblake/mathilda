/* Shared linear-algebra helpers used across the src/linalg/ module and
 * exposed (via linalg.h) to matsol.c, matinv.c, matlstsq.c, mateigen.c.
 *
 * Contains:
 *   get_tensor_dims    -- rank / shape probe for nested-List tensors
 *   flatten_tensor     -- row-major flatten into a pre-allocated array
 *   exact_div_wrapper  -- polynomial exact-division wrapper used by
 *                         RowReduce / LinearSolve to keep intermediate
 *                         entries simplified.
 */

#include "linalg.h"
#include "eval.h"
#include "expand.h"
#include "poly.h"
#include "linsolve.h"   /* matsol_is_inexact -- a zero keeps its numerator's exactness */
#include "sym_names.h"
#include <stdlib.h>

int get_tensor_dims(Expr* e, int64_t* dims) {
    if (e->type != EXPR_FUNCTION || e->data.function.head->type != EXPR_SYMBOL || e->data.function.head->data.symbol.name != SYM_List) {
        return 0; // rank 0
    }
    int64_t len = e->data.function.arg_count;
    dims[0] = len;
    if (len == 0) return 1;

    int sub_rank = get_tensor_dims(e->data.function.args[0], dims + 1);
    for (int64_t i = 1; i < len; i++) {
        int64_t cur_dims[64];
        int cur_rank = get_tensor_dims(e->data.function.args[i], cur_dims);
        if (cur_rank != sub_rank) return -1; // jagged
        for (int j = 0; j < sub_rank; j++) {
            if (cur_dims[j] != dims[j + 1]) return -1; // jagged
        }
    }
    return sub_rank + 1;
}

void flatten_tensor(Expr* e, Expr** flat, size_t* idx) {
    if (e->type == EXPR_FUNCTION && e->data.function.head->type == EXPR_SYMBOL && e->data.function.head->data.symbol.name == SYM_List) {
        for (size_t i = 0; i < e->data.function.arg_count; i++) {
            flatten_tensor(e->data.function.args[i], flat, idx);
        }
    } else {
        flat[(*idx)++] = expr_copy(e);
    }
}

/* Exact division over radical generators.  collect_variables enrols every
 * power of a radical's base as its own generator -- c, Sqrt[c] and c^(3/2) are
 * three independent variables to exact_poly_div -- so a quotient that exists in
 * Q[a, b, Sqrt[c]] (5440 a^3 c^2 / (-10 Sqrt[c]) = -544 a^3 c^(3/2)) is not
 * found, and Bareiss elimination over radical entries then carried an
 * unreduced num/den into its pivot tests: RowReduce, MatrixRank, NullSpace and
 * LinearSolve over Q(params)(Sqrt[c]) reported a consistent system as
 * inconsistent.  Rewrite each radical base B with exponent denominator m as
 * g^m (poly_find_radical_gen, the substitution Together and Cancel use), so
 * the ring is a free polynomial ring in g, divide there, and map back.  The
 * back map g -> B^(1/m) is a ring homomorphism, so a quotient found in the free
 * ring is a quotient of the original; NULL when there is no radical base or
 * the division is not exact there either. */
/* AlgebraicNumber atoms are not free generators (a product of two collapses to
 * one), so the free-ring argument below does not apply to them; their matrices
 * take the number-field RREF in any case. */
static bool has_algebraic_number(const Expr* e) {
    if (e->type != EXPR_FUNCTION) return false;
    if (e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_AlgebraicNumber) return true;
    for (size_t i = 0; i < e->data.function.arg_count; i++)
        if (has_algebraic_number(e->data.function.args[i])) return true;
    return false;
}

static bool radical_base_is_constant(Expr* base) {
    size_t vc = 0, vcap = 4;
    Expr** vs = malloc(sizeof(Expr*) * vcap);
    collect_variables(base, &vs, &vc, &vcap);
    for (size_t i = 0; i < vc; i++) expr_free(vs[i]);
    free(vs);
    return vc == 0;
}

static Expr* exact_div_radical_gens(Expr* exp_num, Expr* exp_den) {
    enum { MAXG = 8 };
    Expr* bases[MAXG]; Expr* atoms[MAXG]; int64_t ms[MAXG]; char* gens[MAXG];
    int ng = 0;
    Expr* pair = expr_new_function(expr_new_symbol(SYM_List),
                                   (Expr*[]){expr_copy(exp_num), expr_copy(exp_den)}, 2);
    while (ng < MAXG) {
        Expr* base = NULL; Expr* atom = NULL; int64_t m = 1;
        if (!poly_find_radical_gen(pair, &base, &atom, &m)) break;
        /* A constant base (Sqrt[2]) is not free: g = 2^(1/2) still satisfies
         * g^2 = 2, which the free ring forgets, so a quotient that exists in
         * Q(Sqrt[2])[x] is found for some entries and not others, and the
         * elimination mixes reduced and unreduced pivots.  Leave such systems
         * to the plain division. */
        if (radical_base_is_constant(base)) {
            expr_free(base); expr_free(atom);
            for (int k = 0; k < ng; k++) { expr_free(bases[k]); expr_free(atoms[k]); free(gens[k]); }
            expr_free(pair);
            return NULL;
        }
        char* gen = poly_make_fresh_gen(pair);
        Expr* sub = poly_subst_radical_to_gen(pair, base, atom, m, gen);
        expr_free(pair);
        pair = expr_expand(sub);
        expr_free(sub);
        bases[ng] = base; atoms[ng] = atom; ms[ng] = m; gens[ng] = gen;
        ng++;
    }
    Expr* res = NULL;
    if (ng > 0 && pair->type == EXPR_FUNCTION && pair->data.function.arg_count == 2) {
        Expr* n2 = pair->data.function.args[0];
        Expr* d2 = pair->data.function.args[1];
        size_t v_count = 0, v_cap = 16;
        Expr** vars = malloc(sizeof(Expr*) * v_cap);
        collect_variables(n2, &vars, &v_count, &v_cap);
        collect_variables(d2, &vars, &v_count, &v_cap);
        if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);
        Expr* q = exact_poly_div(n2, d2, vars, v_count);
        for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
        free(vars);
        if (q) {
            for (int k = ng - 1; k >= 0; k--) {
                Expr* back = poly_subst_radical_from_gen(q, bases[k], atoms[k], ms[k], gens[k]);
                expr_free(q);
                q = back;
            }
            res = expr_expand(q);
            expr_free(q);
            /* self-check in the original ring: q den - num must vanish */
            Expr* chk = eval_and_free(expr_new_function(expr_new_symbol(SYM_Plus), (Expr*[]){
                eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){expr_copy(res), expr_copy(exp_den)}, 2)),
                eval_and_free(expr_new_function(expr_new_symbol(SYM_Times),
                    (Expr*[]){expr_new_integer(-1), expr_copy(exp_num)}, 2))}, 2));
            Expr* chk_e = expr_expand(chk);
            expr_free(chk);
            if (!is_zero_poly(chk_e)) { expr_free(res); res = NULL; }
            expr_free(chk_e);
        }
    }
    for (int k = 0; k < ng; k++) { expr_free(bases[k]); expr_free(atoms[k]); free(gens[k]); }
    expr_free(pair);
    return res;
}

Expr* exact_div_wrapper(Expr* num, Expr* den) {
    /* A zero numerator divides to zero -- but at the NUMERATOR's own exactness.
     * Returning the exact Integer 0 unconditionally put exact zeros into Real
     * results: Inverse[{{2., 0.}, {0., 4.}}] answered {{0.5, 0}, {0, 0.25}}
     * where Mathematica gives all Reals (0./5 is 0., not 0). The heads are the
     * visible half of it; the other half is that a matrix of mixed exact and
     * inexact entries cannot be packed, so every operation downstream of an
     * inverse or a row reduction with a structural zero in it falls off the
     * buffer path. */
    if (is_zero_poly(num))
        return matsol_is_inexact(num) ? expr_new_real(0.0) : expr_new_integer(0);
    if (den->type == EXPR_INTEGER && den->data.integer == 1) return expr_expand(num);

    Expr* exp_num = expr_expand(num);
    Expr* exp_den = expr_expand(den);

    size_t v_count = 0, v_cap = 16;
    Expr** vars = malloc(sizeof(Expr*) * v_cap);
    collect_variables(exp_num, &vars, &v_count, &v_cap);
    collect_variables(exp_den, &vars, &v_count, &v_cap);
    if (v_count > 0) qsort(vars, v_count, sizeof(Expr*), compare_expr_ptrs);

    /* Two constants, one of them inexact: plain division keeps the inexactness.
     * Exact polynomial division in no variables divides structurally, so
     * 1./1. came back as the exact 1 -- RowReduce[{{1., 2., 1.}, {2., 4., 2.}}]
     * answered {{1.0, 2.0, 1}, ...}, a result of two heads from an input of
     * one. Falls through to the Times/Power branch below. */
    Expr* res = (v_count == 0 && (matsol_is_inexact(num) || matsol_is_inexact(den)))
                    ? NULL
                    : exact_poly_div(exp_num, exp_den, vars, v_count);

    for (size_t i = 0; i < v_count; i++) expr_free(vars[i]);
    free(vars);

    if (!res && v_count > 0 && !matsol_is_inexact(num) && !matsol_is_inexact(den)
        && !has_algebraic_number(exp_num) && !has_algebraic_number(exp_den))
        res = exact_div_radical_gens(exp_num, exp_den);

    if (res) {
        expr_free(exp_num);
        expr_free(exp_den);
        return res;
    } else {
        Expr* t = eval_and_free(expr_new_function(expr_new_symbol(SYM_Power), (Expr*[]){exp_den, expr_new_integer(-1)}, 2));
        Expr* r = eval_and_free(expr_new_function(expr_new_symbol(SYM_Times), (Expr*[]){exp_num, t}, 2));
        return r;
    }
}
