/*
 * solvelinsys.c
 *
 * Linear-system specialist for `Solve` (src/solve.c).  Also reachable
 * directly as the context-qualified builtin `Solve`SolveLinearSystem`.
 *
 * Algorithm
 * ---------
 * Each equation `lhs_i == rhs_i` is canonicalised to `poly_i := lhs_i
 * - rhs_i`.  We then assert `poly_i` is *affine* (degree at most 1) in
 * every input variable -- otherwise the specialist refuses and the
 * caller leaves Solve unevaluated.  Affinity check:
 *
 *   1. is_polynomial(poly_i, vars)            -- structural test
 *   2. coeff_j = get_coeff(poly_i, vars[j], 1) for each j;
 *      verify each coeff_j is free of every var (no cross terms).
 *   3. residual = poly_i - sum_j coeff_j * vars[j];
 *      verify residual is free of every var (no higher-degree terms).
 *
 * Coefficients populate an m x (n+1) augmented matrix M with the
 * variable columns *reversed* -- M[i][0] is the coefficient of
 * vars[n-1], M[i][n-1] of vars[0], and M[i][n] = -residual (so the
 * row represents `sum_j coeff_j * v_j = -residual`).
 *
 * The reversed-column trick is what produces Mathematica's
 * `Solve::svars` convention: with vars = {x, y} and one equation,
 * standard left-to-right Gaussian elimination naturally pivots on
 * the reversed-first column = original `y`, leaving `x` free.
 *
 * We then run Gauss--Jordan to reduced row echelon form, with one
 * subtlety for symbolic pivot selection: we prefer concretely
 * non-zero pivots (Integer / Rational / Real) over symbolic ones,
 * but accept any expression that does not simplify to zero via
 * Cancel[Together[.]].  After reduction:
 *
 *   - rank == #vars   ->  unique solution; read off rules.
 *   - rank <  #vars   ->  free variables exist; emit Solve::svars;
 *                         each pivot variable becomes a rule whose RHS
 *                         is the pivot row's augmented column minus
 *                         the free-column terms re-expressed in the
 *                         free variables.
 *   - any zero row    ->  if the augmented column is non-zero in that
 *                         row, the system is inconsistent -> return {}.
 *
 * Domain filter (applied after solving):
 *   - Integers : every emitted rule's RHS must be a concrete Integer /
 *                BigInt; otherwise the whole solution is dropped.
 *   - Reals    : drop a solution if any RHS contains a Complex[_, _]
 *                head (a structurally-non-real value).
 *   - Complexes: default; no filter.
 *
 * Memory: every entry of the augmented matrix is freshly owned; the
 * helpers below maintain that invariant.  All paths through this file
 * must either free or return ownership of every Expr* they touch.
 */

#include "solvelinsys.h"

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "arithmetic.h"
#include "attr.h"
#include "eval.h"
#include "expr.h"
#include "internal.h"
#include "poly.h"
#include "sym_names.h"
#include "symtab.h"

/* ------------------------------------------------------------------ *
 *  Expression-building shorthand (same conventions as solvepoly.c).   *
 *  The mk_* helpers take ownership of their pointer arguments.        *
 * ------------------------------------------------------------------ */

static Expr* mk_int(int64_t v) { return expr_new_integer(v); }
static Expr* mk_sym(const char* name) { return expr_new_symbol(name); }

static Expr* mk_fn2(const char* head, Expr* a, Expr* b) {
    return expr_new_function(mk_sym(head), (Expr*[]){ a, b }, 2);
}
static Expr* mk_pow(Expr* base, Expr* exp) { return mk_fn2("Power", base, exp); }
static Expr* mk_neg(Expr* e) { return mk_fn2("Times", mk_int(-1), e); }
static Expr* mk_rule(Expr* lhs, Expr* rhs) { return mk_fn2("Rule", lhs, rhs); }

static Expr* mk_list(Expr** args, size_t n) {
    return expr_new_function(mk_sym("List"), args, n);
}

/* ------------------------------------------------------------------ *
 *  Zero-testing and entry simplification.                             *
 * ------------------------------------------------------------------ */

/* True iff `e` is structurally the integer 0. */
static bool is_int_zero(const Expr* e) {
    return e && e->type == EXPR_INTEGER && e->data.integer == 0;
}

/* Normalise `e` to a single rational/algebraic form via
 * Cancel[Together[e]].  Takes ownership of `e`; returns a freshly
 * owned canonical form.  This is what makes the symbolic zero-test
 * below decisive for the cases that matter for elimination
 * (`a*x - a*x`, `(p/q) - p*Power[q,-1]`, ...). */
static Expr* simplify_entry(Expr* e) {
    Expr* tog = internal_together((Expr*[]){ e }, 1);
    Expr* can = internal_cancel((Expr*[]){ tog }, 1);
    return can;
}

/* True iff `e` simplifies (via Cancel[Together[.]]) to the integer 0.
 * Borrowed argument. */
static bool simplifies_to_zero(const Expr* e) {
    if (!e) return true;
    if (is_int_zero(e)) return true;
    if (e->type == EXPR_INTEGER) return false;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    /* For everything else, run the full simplifier. */
    Expr* tmp = simplify_entry(expr_copy((Expr*)e));
    bool zero = is_int_zero(tmp);
    expr_free(tmp);
    return zero;
}

/* Concrete-numeric pivot detector: true for an entry where we can be
 * certain it is non-zero by inspecting its head alone (no need for
 * symbolic case-splitting).  Used to prefer Integer/Rational/Real
 * pivots over symbolic ones during pivot selection. */
static bool is_concrete_nonzero(const Expr* e) {
    if (!e) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer != 0;
    if (e->type == EXPR_REAL)    return e->data.real != 0.0;
    if (e->type == EXPR_FUNCTION
        && e->data.function.head->type == EXPR_SYMBOL
        && e->data.function.head->data.symbol.name == SYM_Rational) {
        /* Mathilda canonicalises Rational[0, _] -> 0, so any
         * Rational[_, _] surviving here is non-zero. */
        return true;
    }
    /* BigInt is reachable via SYM lookups; treat conservatively. */
    return false;
}

/* True iff `v` is provably a concrete integer.  Used by the
 * `Integers` domain filter. */
static bool is_concrete_integer(const Expr* v) {
    if (!v) return false;
    return v->type == EXPR_INTEGER;
}

/* True iff `v` syntactically contains a `Complex[_, _]` node anywhere.
 * Heuristic used by the `Reals` filter to drop provably non-real
 * solutions. */
static bool contains_complex_head(const Expr* v) {
    if (!v) return false;
    if (v->type != EXPR_FUNCTION) return false;
    if (v->data.function.head->type == EXPR_SYMBOL
        && v->data.function.head->data.symbol.name == SYM_Complex) {
        return true;
    }
    for (size_t i = 0; i < v->data.function.arg_count; i++) {
        if (contains_complex_head(v->data.function.args[i])) return true;
    }
    return false;
}

/* ------------------------------------------------------------------ *
 *  Affine-extraction: turn one equation into a coefficient row.        *
 * ------------------------------------------------------------------ */

/* Given a polynomial-form expression `poly` (= lhs - rhs of an
 * equation) and a list of `n` variable symbols `vars[]`, extract the
 * `n` linear coefficients and the constant residual.  Returns true on
 * success with `coeffs_out[j]` (freshly owned) set for each j and
 * `*constant_out` (freshly owned) holding the var-free residual.
 *
 * Returns false (with everything cleaned up) when `poly` is non-affine
 * in the vars (cross terms, higher degrees, transcendental dependence).
 *
 * Linearity is verified two ways:
 *   1. Each coeff_j must be free of every var (rules out cross terms
 *      such as `a*x*y` whose `Coefficient[..., x, 1]` would be `a*y`).
 *   2. residual = poly - sum_j coeff_j * v_j must be free of every var
 *      (rules out higher-degree terms such as `x^2`). */
static bool extract_affine(Expr* poly, Expr** vars, size_t n,
                           Expr** coeffs_out, Expr** constant_out) {
    if (!is_polynomial(poly, vars, n)) return false;

    for (size_t j = 0; j < n; j++) coeffs_out[j] = NULL;

    for (size_t j = 0; j < n; j++) {
        Expr* c = get_coeff(poly, vars[j], 1);
        if (!c) {
            for (size_t k = 0; k < j; k++) expr_free(coeffs_out[k]);
            return false;
        }
        for (size_t k = 0; k < n; k++) {
            if (contains_any_symbol_from(c, vars[k])) {
                expr_free(c);
                for (size_t kk = 0; kk < j; kk++) expr_free(coeffs_out[kk]);
                return false;
            }
        }
        coeffs_out[j] = c;
    }

    /* residual = poly - sum_j c_j * v_j.  Build as Plus[poly, -c0*v0,
     * -c1*v1, ...] and evaluate-then-simplify.  */
    Expr** terms = (Expr**)malloc(sizeof(Expr*) * (n + 1));
    terms[0] = expr_copy(poly);
    for (size_t j = 0; j < n; j++) {
        Expr* prod = mk_fn2("Times", expr_copy(coeffs_out[j]),
                                     expr_copy(vars[j]));
        terms[j + 1] = mk_neg(prod);
    }
    Expr* sum = expr_new_function(mk_sym("Plus"), terms, n + 1);
    free(terms);
    Expr* residual = simplify_entry(eval_and_free(sum));

    for (size_t k = 0; k < n; k++) {
        if (contains_any_symbol_from(residual, vars[k])) {
            expr_free(residual);
            for (size_t kk = 0; kk < n; kk++) expr_free(coeffs_out[kk]);
            return false;
        }
    }

    *constant_out = residual;
    return true;
}

/* ------------------------------------------------------------------ *
 *  Augmented matrix and Gauss--Jordan elimination.                    *
 * ------------------------------------------------------------------ */

/* An m x (n+1) augmented matrix where each cell owns its Expr*.
 * Columns 0..n-1 are coefficient columns (vars in REVERSED order:
 * M[i][0] is the coefficient of vars[n-1]).  Column n is the negated
 * constant (`sum_j c_j * v_j = -const`).  */
typedef struct {
    Expr*** rows;
    size_t m, n;
} AugMat;

static void am_init(AugMat* A, size_t m, size_t n) {
    A->m = m; A->n = n;
    A->rows = (Expr***)malloc(sizeof(Expr**) * (m ? m : 1));
    for (size_t i = 0; i < m; i++) {
        A->rows[i] = (Expr**)malloc(sizeof(Expr*) * (n + 1));
        for (size_t j = 0; j <= n; j++) A->rows[i][j] = NULL;
    }
}

static void am_free(AugMat* A) {
    if (!A->rows) return;
    for (size_t i = 0; i < A->m; i++) {
        for (size_t j = 0; j <= A->n; j++) {
            if (A->rows[i][j]) expr_free(A->rows[i][j]);
        }
        free(A->rows[i]);
    }
    free(A->rows);
    A->rows = NULL;
}

static void am_swap_rows(AugMat* A, size_t a, size_t b) {
    if (a == b) return;
    Expr** tmp = A->rows[a];
    A->rows[a] = A->rows[b];
    A->rows[b] = tmp;
}

/* Replace A[i][j] with simplify(new_val).  Takes ownership of new_val. */
static void am_set(AugMat* A, size_t i, size_t j, Expr* new_val) {
    if (A->rows[i][j]) expr_free(A->rows[i][j]);
    A->rows[i][j] = simplify_entry(new_val);
}

/* Run Gauss--Jordan elimination to reduced row echelon form.  Writes
 * the pivot row for each column into pivot_row_for_col[]; columns
 * without a pivot (= free variables) receive SIZE_MAX.
 *
 * Returns -1 if any inconsistent row of the form `0 = c` (c ≠ 0) is
 * detected, else the rank (number of pivots).
 *
 * Pivot selection: among rows [r0, m) for the current pivot column,
 * prefer a concretely-nonzero pivot (Integer / Rational / Real) for
 * numerical robustness; otherwise accept the first entry that does
 * not simplify to zero.  All-zero column => free variable. */
static int gaussian_reduce(AugMat* A, size_t* pivot_row_for_col) {
    for (size_t j = 0; j < A->n; j++) pivot_row_for_col[j] = SIZE_MAX;

    size_t pivot_row = 0;
    for (size_t j = 0; j < A->n && pivot_row < A->m; j++) {
        /* Pivot selection. */
        size_t pick = SIZE_MAX;
        size_t fallback = SIZE_MAX;
        for (size_t r = pivot_row; r < A->m; r++) {
            if (simplifies_to_zero(A->rows[r][j])) continue;
            if (is_concrete_nonzero(A->rows[r][j])) { pick = r; break; }
            if (fallback == SIZE_MAX) fallback = r;
        }
        if (pick == SIZE_MAX) pick = fallback;
        if (pick == SIZE_MAX) continue;            /* free variable column */

        am_swap_rows(A, pivot_row, pick);
        pivot_row_for_col[j] = pivot_row;

        /* Normalise pivot row so A[pivot_row][j] == 1. */
        Expr* inv_pivot = eval_and_free(
            mk_pow(expr_copy(A->rows[pivot_row][j]), mk_int(-1)));
        for (size_t k = j; k <= A->n; k++) {
            Expr* prod = eval_and_free(mk_fn2("Times",
                expr_copy(A->rows[pivot_row][k]),
                expr_copy(inv_pivot)));
            am_set(A, pivot_row, k, prod);
        }
        expr_free(inv_pivot);

        /* Eliminate column j in every other row. */
        for (size_t r = 0; r < A->m; r++) {
            if (r == pivot_row) continue;
            if (simplifies_to_zero(A->rows[r][j])) continue;
            Expr* factor = expr_copy(A->rows[r][j]);
            for (size_t k = j; k <= A->n; k++) {
                Expr* delta = eval_and_free(mk_fn2("Times",
                    expr_copy(factor),
                    expr_copy(A->rows[pivot_row][k])));
                Expr* new_val = eval_and_free(mk_fn2("Plus",
                    expr_copy(A->rows[r][k]),
                    mk_neg(delta)));
                am_set(A, r, k, new_val);
            }
            expr_free(factor);
        }
        pivot_row++;
    }

    /* Inconsistency check: any row past the last pivot row with a
     * non-zero augmented column => `0 = c` with c ≠ 0. */
    int rank = (int)pivot_row;
    for (size_t r = (size_t)rank; r < A->m; r++) {
        if (!simplifies_to_zero(A->rows[r][A->n])) {
            return -1;
        }
    }
    return rank;
}

/* ------------------------------------------------------------------ *
 *  svars warning.                                                     *
 * ------------------------------------------------------------------ */

/* Emit `Solve::svars` to stderr.  Throttled to once per distinct
 * input hash, mirroring the `Solve::optx` idiom in solve.c. */
static void warn_svars(uint64_t input_hash) {
    /* Internal decision procedures (e.g. the Risch structure-theorem Q-span
     * solve) legitimately hand SolveAlways underdetermined systems; honour the
     * shared warnings-mute counter so those probes stay quiet. */
    if (arith_warnings_muted()) return;
    static uint64_t last_warned_hash = 0;
    if (input_hash == last_warned_hash) return;
    last_warned_hash = input_hash;
    fprintf(stderr,
        "Solve::svars: Equations may not give solutions for all "
        "\"solve\" variables.\n");
}

/* ------------------------------------------------------------------ *
 *  Public entry.                                                      *
 * ------------------------------------------------------------------ */

Expr* solvelinsys_solve_linear_system(Expr* equations,
                                      Expr* vars,
                                      Expr* dom) {
    if (!equations || !vars) return NULL;

    /* `vars` must be a List of symbols (>= 1). */
    if (vars->type != EXPR_FUNCTION
        || vars->data.function.head->type != EXPR_SYMBOL
        || vars->data.function.head->data.symbol.name != SYM_List) {
        return NULL;
    }
    size_t n = vars->data.function.arg_count;
    if (n == 0) return NULL;
    Expr** var_arr = vars->data.function.args;
    for (size_t j = 0; j < n; j++) {
        if (var_arr[j]->type != EXPR_SYMBOL) return NULL;
    }
    /* Reject duplicate variables -- semantics would be ambiguous and
     * Mathematica also refuses. */
    for (size_t a = 0; a < n; a++) {
        for (size_t b = a + 1; b < n; b++) {
            if (var_arr[a]->data.symbol.name == var_arr[b]->data.symbol.name) {
                return NULL;
            }
        }
    }

    /* `equations` is one of: Equal[_, _], And[Equal, ...], List[Equal, ...]. */
    Expr* single_holder[1];
    Expr** eq_arr;
    size_t m;
    if (equations->type == EXPR_FUNCTION
        && equations->data.function.head->type == EXPR_SYMBOL
        && (equations->data.function.head->data.symbol.name == SYM_List
            || equations->data.function.head->data.symbol.name == SYM_And)) {
        eq_arr = equations->data.function.args;
        m = equations->data.function.arg_count;
    } else if (equations->type == EXPR_FUNCTION
        && equations->data.function.head->type == EXPR_SYMBOL
        && equations->data.function.head->data.symbol.name == SYM_Equal) {
        single_holder[0] = equations;
        eq_arr = single_holder;
        m = 1;
    } else {
        return NULL;
    }

    /* Empty system -- vacuously true, every assignment satisfies.
     * Mirror Mathematica: Solve[True, vars] returns { {} }. */
    if (m == 0) {
        Expr* empty = mk_list(NULL, 0);
        return mk_list((Expr*[]){ empty }, 1);
    }

    /* Domain. */
    bool reals_only = false;
    bool integers_only = false;
    if (dom) {
        if (dom->type != EXPR_SYMBOL) return NULL;
        const char* d = dom->data.symbol.name;
        if (d == SYM_Reals) reals_only = true;
        else if (d == SYM_Integers) {
            reals_only = true;
            integers_only = true;
        } else if (d != SYM_Complexes) {
            return NULL;
        }
    }

    /* Build the augmented matrix, one row per equation. */
    AugMat A;
    am_init(&A, m, n);

    for (size_t i = 0; i < m; i++) {
        Expr* eq = eq_arr[i];
        if (eq->type != EXPR_FUNCTION
            || eq->data.function.head->type != EXPR_SYMBOL
            || eq->data.function.head->data.symbol.name != SYM_Equal
            || eq->data.function.arg_count != 2) {
            am_free(&A);
            return NULL;
        }
        Expr* lhs = eq->data.function.args[0];
        Expr* rhs = eq->data.function.args[1];

        /* poly = lhs - rhs, expanded to standard form so coefficients
         * surface cleanly (e.g. `a (x + y) + b` becomes `a x + a y + b`). */
        Expr* diff = eval_and_free(mk_fn2("Plus",
            expr_copy(lhs), mk_neg(expr_copy(rhs))));
        Expr* poly = internal_expand((Expr*[]){ diff }, 1);

        Expr** coeffs = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
        Expr* constant = NULL;
        bool ok = extract_affine(poly, var_arr, n, coeffs, &constant);
        expr_free(poly);
        if (!ok) {
            free(coeffs);
            am_free(&A);
            return NULL;
        }

        /* Row i, reversed-column layout: M[i][k] = coeff of vars[n-1-k].
         * Augmented entry M[i][n] = -constant. */
        for (size_t k = 0; k < n; k++) {
            A.rows[i][k] = simplify_entry(coeffs[n - 1 - k]);
        }
        free(coeffs);
        A.rows[i][n] = simplify_entry(mk_neg(constant));
    }

    /* Reduce to RREF. */
    size_t* pivot_row_for_col = (size_t*)malloc(sizeof(size_t) * n);
    int rank = gaussian_reduce(&A, pivot_row_for_col);

    if (rank < 0) {
        free(pivot_row_for_col);
        am_free(&A);
        return mk_list(NULL, 0);    /* inconsistent  ->  {} */
    }

    /* Build solution rules in original variable order.  A column k
     * (reversed) corresponds to original variable index (n - 1 - k). */
    bool* var_is_pivot = (bool*)calloc(n, sizeof(bool));
    for (size_t k = 0; k < n; k++) {
        if (pivot_row_for_col[k] != SIZE_MAX) {
            var_is_pivot[n - 1 - k] = true;
        }
    }
    bool any_free = false;
    for (size_t j = 0; j < n; j++) if (!var_is_pivot[j]) { any_free = true; break; }

    Expr** rules = (Expr**)malloc(sizeof(Expr*) * (n ? n : 1));
    size_t rules_count = 0;

    for (size_t j = 0; j < n; j++) {
        if (!var_is_pivot[j]) continue;
        size_t k = n - 1 - j;                  /* reversed column index */
        size_t pr = pivot_row_for_col[k];

        /* val = A[pr][n] - sum_{l != k, l is free} A[pr][l] * vars[n-1-l].
         * (Pivot-pivot off-diagonals are zero in RREF, so we only need
         * to sweep free columns -- but iterating across all l != k and
         * skipping zeros is just as safe and is robust to any RREF
         * shortcuts we may make later.) */
        Expr* val = expr_copy(A.rows[pr][A.n]);
        for (size_t l = 0; l < A.n; l++) {
            if (l == k) continue;
            if (simplifies_to_zero(A.rows[pr][l])) continue;
            Expr* term = eval_and_free(mk_fn2("Times",
                expr_copy(A.rows[pr][l]),
                expr_copy(var_arr[n - 1 - l])));
            val = eval_and_free(mk_fn2("Plus", val, mk_neg(term)));
        }
        val = simplify_entry(val);

        /* Domain filter: reject the whole solution if this value
         * violates the requested domain. */
        bool reject = false;
        if (integers_only && !is_concrete_integer(val)) reject = true;
        else if (reals_only && contains_complex_head(val)) reject = true;
        if (reject) {
            expr_free(val);
            for (size_t r = 0; r < rules_count; r++) expr_free(rules[r]);
            free(rules);
            free(var_is_pivot);
            free(pivot_row_for_col);
            am_free(&A);
            return mk_list(NULL, 0);
        }

        rules[rules_count++] = mk_rule(expr_copy(var_arr[j]), val);
    }

    am_free(&A);
    free(pivot_row_for_col);
    free(var_is_pivot);

    if (any_free) {
        /* Use the input's hash so repeated identical calls don't spam. */
        Expr* probe = mk_fn2("SolveLinearSystem",
            expr_copy(equations), expr_copy(vars));
        warn_svars(expr_hash(probe));
        expr_free(probe);
    }

    /* expr_new_function memcpys the args[] array (src/expr.c:88-93),
     * so the heap-allocated `rules` buffer is no longer needed once
     * mk_list has copied its pointers. */
    Expr* rule_list = mk_list(rules, rules_count);
    free(rules);
    return mk_list((Expr*[]){ rule_list }, 1);
}

/* ------------------------------------------------------------------ *
 *  Builtin entry & init.                                              *
 * ------------------------------------------------------------------ */

Expr* builtin_solve_linear_system(Expr* res) {
    if (!res || res->type != EXPR_FUNCTION) return NULL;
    size_t argc = res->data.function.arg_count;
    if (argc < 2 || argc > 3) return NULL;
    Expr* equations = res->data.function.args[0];
    Expr* vars      = res->data.function.args[1];
    Expr* dom       = (argc >= 3) ? res->data.function.args[2] : NULL;
    return solvelinsys_solve_linear_system(equations, vars, dom);
}

void solvelinsys_init(void) {
    symtab_add_builtin("Solve`SolveLinearSystem",
                       builtin_solve_linear_system);
    SymbolDef* def = symtab_get_def("Solve`SolveLinearSystem");
    if (def) def->attributes |= ATTR_PROTECTED;
    symtab_set_docstring("Solve`SolveLinearSystem",
        "Solve`SolveLinearSystem[eqns, vars]\n"
        "Solve`SolveLinearSystem[eqns, vars, dom]\n"
        "\tThe linear-system specialist used by Solve.  Solves the\n"
        "\tconjunction (And, List, or single Equal) `eqns` for the\n"
        "\tvariables `vars` via Gauss--Jordan elimination.\n"
        "\n"
        "\tEach equation lhs == rhs is canonicalised to lhs - rhs and\n"
        "\tmust be affine in `vars` -- otherwise the call is left\n"
        "\tunevaluated.\n"
        "\n"
        "\tReturns:\n"
        "\t  unique solution        {{v1 -> e1, v2 -> e2, ...}}\n"
        "\t  inconsistent system    {}\n"
        "\t  underdetermined system {{pivot_vars -> exprs}}\n"
        "\t                         (emits Solve::svars; free variables\n"
        "\t                         produce no rule)\n"
        "\n"
        "\tdom = Complexes (default), Reals, or Integers.  Reals drops\n"
        "\tsolutions whose RHS contains a non-real Complex[_, _]; Integers\n"
        "\tdrops any solution whose RHS is not a provably concrete integer.");
}
