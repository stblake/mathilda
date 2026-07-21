#include "list_common.h"
#include "matrixq.h"
#include "ndarray.h"

/* An NDArray is a distinct atomic node, so the List-based matrix predicates
 * below would otherwise answer False for a genuine packed matrix. Route any
 * NDArray argument through the exact List path (ndarray_delist_and_reeval
 * rebuilds the call with NDArrays converted to nested Lists and re-evaluates —
 * no recursion, since the rebuilt call has no NDArray args). */
#define ND_MATRIXQ_DELIST(res)                                                  \
    do {                                                                        \
        if ((res)->data.function.arg_count >= 1 &&                              \
            is_ndarray((res)->data.function.args[0]))                           \
            return ndarray_delist_and_reeval(res);                             \
    } while (0)

/* --- HermitianMatrixQ ----------------------------------------------------
 *
 * A matrix m is Hermitian iff m == ConjugateTranspose[m], equivalently
 * m[i,j] == Conjugate[m[j,i]] for every pair (i,j).  Diagonal entries
 * must satisfy m[i,i] == Conjugate[m[i,i]], i.e. be self-conjugate
 * (purely real for numeric inputs).
 *
 * The default test is "explicit" (structural) matching: we accept a pair
 * (a, b) as conjugate-mirrored when any of the following holds --
 *
 *   (1) a is Conjugate[c] with c structurally equal to b;
 *   (2) b is Conjugate[c] with c structurally equal to a;
 *   (3) evaluating Conjugate[b] yields a structurally.
 *
 * Branches (1) and (2) catch the symbolic patterns (Conjugate[a], a) and
 * (a, Conjugate[a]) because our Conjugate builtin does NOT fold
 * Conjugate[Conjugate[x]] back to x for symbolic x.  Branch (3) covers
 * the fully numeric case where Conjugate evaluates to a concrete value.
 *
 * Options:
 *   - SameTest -> f : pairs (a, b) are accepted iff f[a, Conjugate[b]]
 *     evaluates to True.  Conjugate[b] is computed once per pair so the
 *     user-supplied predicate sees the actual conjugated entry, mirroring
 *     Mathematica's documented behaviour.
 *   - Tolerance -> t : pairs (a, b) are accepted iff
 *     Abs[a - Conjugate[b]] <= t evaluates to True.
 *
 * SameTest and Tolerance defaults of Automatic fall through to the
 * structural test, which is correct for both symbolic and exact-numeric
 * matrices and degrades to bit-identical comparison for machine reals.
 */

/* True when `n` is an evaluated Conjugate[x] node. */
static bool is_conjugate_node(Expr* n) {
    return n->type == EXPR_FUNCTION &&
           n->data.function.head->type == EXPR_SYMBOL &&
           n->data.function.head->data.symbol.name == SYM_Conjugate &&
           n->data.function.arg_count == 1;
}

static Expr* eval_conjugate_of(Expr* e) {
    Expr** args = malloc(sizeof(Expr*) * 1);
    args[0] = expr_copy(e);
    Expr* call = expr_new_function(expr_new_symbol(SYM_Conjugate), args, 1);
    free(args);
    return eval_and_free(call);
}

static bool hermitian_pair_structural(Expr* a, Expr* b) {
    /* (1) a == Conjugate[c] structurally and c == b. */
    if (is_conjugate_node(a) && expr_eq(a->data.function.args[0], b)) return true;
    /* (2) b == Conjugate[c] structurally and c == a. */
    if (is_conjugate_node(b) && expr_eq(b->data.function.args[0], a)) return true;
    /* (3) Evaluate Conjugate[b] and compare structurally. */
    Expr* cb = eval_conjugate_of(b);
    bool eq = expr_eq(a, cb);
    expr_free(cb);
    return eq;
}

static bool hermitian_pair_sametest(Expr* a, Expr* b, Expr* test) {
    /* Build test[a, Conjugate[b]] and check for True. */
    Expr* cb = eval_conjugate_of(b);
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(a);
    args[1] = cb;  /* ownership transferred into the call below */
    Expr* call = expr_new_function(expr_copy(test), args, 2);
    free(args);
    Expr* result = eval_and_free(call);
    bool ok = (result->type == EXPR_SYMBOL &&
               result->data.symbol.name == SYM_True);
    expr_free(result);
    return ok;
}

static bool hermitian_pair_tolerance(Expr* a, Expr* b, Expr* tol) {
    /* Build LessEqual[Abs[a - Conjugate[b]], tol] and check for True. */
    Expr* cb = eval_conjugate_of(b);
    Expr** sub_args = malloc(sizeof(Expr*) * 2);
    sub_args[0] = expr_copy(a);
    sub_args[1] = cb;
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Subtract), sub_args, 2);
    free(sub_args);
    Expr* diff_e = eval_and_free(diff);

    Expr** abs_args = malloc(sizeof(Expr*) * 1);
    abs_args[0] = diff_e;
    Expr* abs_call = expr_new_function(expr_new_symbol(SYM_Abs), abs_args, 1);
    free(abs_args);
    Expr* abs_e = eval_and_free(abs_call);

    Expr** le_args = malloc(sizeof(Expr*) * 2);
    le_args[0] = abs_e;
    le_args[1] = expr_copy(tol);
    Expr* le_call = expr_new_function(expr_new_symbol(SYM_LessEqual), le_args, 2);
    free(le_args);
    Expr* le_e = eval_and_free(le_call);
    bool ok = (le_e->type == EXPR_SYMBOL &&
               le_e->data.symbol.name == SYM_True);
    expr_free(le_e);
    return ok;
}

Expr* builtin_hermitian_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    ND_MATRIXQ_DELIST(res);
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* m = res->data.function.args[0];

    /* Parse options.  Each trailing arg must be a Rule with one of the
     * recognised option names; any unrecognised option causes us to
     * leave the call unevaluated so user-typed errors surface. */
    Expr* same_test = NULL;
    Expr* tolerance = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!(opt->type == EXPR_FUNCTION &&
              opt->data.function.head->type == EXPR_SYMBOL &&
              opt->data.function.head->data.symbol.name == SYM_Rule &&
              opt->data.function.arg_count == 2 &&
              opt->data.function.args[0]->type == EXPR_SYMBOL)) {
            return NULL;
        }
        const char* name = opt->data.function.args[0]->data.symbol.name;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_SameTest) {
            /* Automatic falls through to the structural test. */
            if (!(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic)) {
                same_test = val;
            }
        } else if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            return NULL;
        }
    }

    /* Must be a non-empty square List of Lists with no deeper nesting. */
    if (!is_listq(m)) return expr_new_symbol(SYM_False);
    size_t n = m->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_False);
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);
        if (row->data.function.arg_count != n) return expr_new_symbol(SYM_False);
        for (size_t j = 0; j < n; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol(SYM_False);
            }
        }
    }

    /* Walk the upper triangle (including diagonal) and check each pair
     * against the chosen predicate.  Walking i <= j is sufficient since
     * the pair test is symmetric under (i,j) <-> (j,i): a == Conj[b]
     * iff b == Conj[a]. */
    for (size_t i = 0; i < n; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = i; j < n; j++) {
            Expr* a = row_i[j];                                /* m[i,j] */
            Expr* b = m->data.function.args[j]->data.function.args[i]; /* m[j,i] */
            bool ok;
            if (same_test != NULL) {
                ok = hermitian_pair_sametest(a, b, same_test);
            } else if (tolerance != NULL) {
                ok = hermitian_pair_tolerance(a, b, tolerance);
            } else {
                ok = hermitian_pair_structural(a, b);
            }
            if (!ok) return expr_new_symbol(SYM_False);
        }
    }

    return expr_new_symbol(SYM_True);
}

/* --- SymmetricMatrixQ ----------------------------------------------------
 *
 * A matrix m is symmetric iff m == Transpose[m], i.e. m[i,j] == m[j,i]
 * for every pair (i, j).  This is the complex-symmetric notion -- no
 * conjugation is applied -- so a complex symmetric matrix need not be
 * Hermitian (and a Hermitian complex matrix need not be symmetric).
 *
 * The default test is structural equality via expr_eq.  Options mirror
 * HermitianMatrixQ:
 *   - SameTest -> f : pairs (a, b) accepted iff f[a, b] evaluates to True.
 *   - Tolerance -> t : pairs accepted iff Abs[a - b] <= t evaluates to True.
 *
 * SameTest/Tolerance defaults of Automatic fall through to the structural
 * test, which is correct for both symbolic and exact-numeric matrices and
 * degrades to bit-identical comparison for machine reals.  The diagonal is
 * always trivially symmetric (m[i,i] == m[i,i]) so we walk only the
 * strict upper triangle.
 */

static bool symmetric_pair_sametest(Expr* a, Expr* b, Expr* test) {
    Expr** args = malloc(sizeof(Expr*) * 2);
    args[0] = expr_copy(a);
    args[1] = expr_copy(b);
    Expr* call = expr_new_function(expr_copy(test), args, 2);
    free(args);
    Expr* result = eval_and_free(call);
    bool ok = (result->type == EXPR_SYMBOL &&
               result->data.symbol.name == SYM_True);
    expr_free(result);
    return ok;
}

static bool symmetric_pair_tolerance(Expr* a, Expr* b, Expr* tol) {
    /* Build LessEqual[Abs[a - b], tol] and check for True. */
    Expr** sub_args = malloc(sizeof(Expr*) * 2);
    sub_args[0] = expr_copy(a);
    sub_args[1] = expr_copy(b);
    Expr* diff = expr_new_function(expr_new_symbol(SYM_Subtract), sub_args, 2);
    free(sub_args);
    Expr* diff_e = eval_and_free(diff);

    Expr** abs_args = malloc(sizeof(Expr*) * 1);
    abs_args[0] = diff_e;
    Expr* abs_call = expr_new_function(expr_new_symbol(SYM_Abs), abs_args, 1);
    free(abs_args);
    Expr* abs_e = eval_and_free(abs_call);

    Expr** le_args = malloc(sizeof(Expr*) * 2);
    le_args[0] = abs_e;
    le_args[1] = expr_copy(tol);
    Expr* le_call = expr_new_function(expr_new_symbol(SYM_LessEqual), le_args, 2);
    free(le_args);
    Expr* le_e = eval_and_free(le_call);
    bool ok = (le_e->type == EXPR_SYMBOL &&
               le_e->data.symbol.name == SYM_True);
    expr_free(le_e);
    return ok;
}

Expr* builtin_symmetric_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    ND_MATRIXQ_DELIST(res);
    size_t argc = res->data.function.arg_count;
    if (argc < 1) return NULL;

    Expr* m = res->data.function.args[0];

    /* Parse options.  Each trailing arg must be a Rule with one of the
     * recognised option names; any unrecognised option causes us to
     * leave the call unevaluated so user-typed errors surface. */
    Expr* same_test = NULL;
    Expr* tolerance = NULL;
    for (size_t i = 1; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        if (!(opt->type == EXPR_FUNCTION &&
              opt->data.function.head->type == EXPR_SYMBOL &&
              opt->data.function.head->data.symbol.name == SYM_Rule &&
              opt->data.function.arg_count == 2 &&
              opt->data.function.args[0]->type == EXPR_SYMBOL)) {
            return NULL;
        }
        const char* name = opt->data.function.args[0]->data.symbol.name;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_SameTest) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic)) {
                same_test = val;
            }
        } else if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            return NULL;
        }
    }

    /* Must be a non-empty square List of Lists with no deeper nesting. */
    if (!is_listq(m)) return expr_new_symbol(SYM_False);
    size_t n = m->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_False);
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);
        if (row->data.function.arg_count != n) return expr_new_symbol(SYM_False);
        for (size_t j = 0; j < n; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol(SYM_False);
            }
        }
    }

    /* Walk strict upper triangle (i < j); the diagonal is trivially
     * symmetric under any reasonable equality test that is reflexive. */
    for (size_t i = 0; i < n; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = i + 1; j < n; j++) {
            Expr* a = row_i[j];                                /* m[i,j] */
            Expr* b = m->data.function.args[j]->data.function.args[i]; /* m[j,i] */
            bool ok;
            if (same_test != NULL) {
                ok = symmetric_pair_sametest(a, b, same_test);
            } else if (tolerance != NULL) {
                ok = symmetric_pair_tolerance(a, b, tolerance);
            } else {
                ok = expr_eq(a, b);
            }
            if (!ok) return expr_new_symbol(SYM_False);
        }
    }

    return expr_new_symbol(SYM_True);
}

/* --- SquareMatrixQ -----------------------------------------------------
 *
 * A matrix m is square iff Dimensions[m] == {n, n}, i.e. it has the
 * same number of rows and columns.  Pure shape test -- no element
 * predicate or option is consulted.  Returns False for non-lists,
 * empty lists, lists whose elements are not all lists of the same
 * length, ragged matrices, rectangular matrices, and higher-rank
 * tensors (any entry that is itself a List).  Single-element matrices
 * {{x}} are square for any x (including symbolic x).
 *
 * Exactly one argument is accepted; any other count emits a
 * Mathematica-compatible `SquareMatrixQ::argx` diagnostic to stderr
 * and leaves the call unevaluated, matching the surface behaviour of
 * other wrong-arity builtins (cf. builtin_conjugate).
 */
Expr* builtin_square_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    ND_MATRIXQ_DELIST(res);
    if (res->data.function.arg_count != 1) {
        size_t n = res->data.function.arg_count;
        fprintf(stderr,
                "SquareMatrixQ::argx: SquareMatrixQ called with %zu "
                "argument%s; 1 argument is expected.\n",
                n, n == 1 ? "" : "s");
        return NULL;
    }

    Expr* m = res->data.function.args[0];

    if (!is_listq(m)) return expr_new_symbol(SYM_False);
    size_t n = m->data.function.arg_count;
    if (n == 0) return expr_new_symbol(SYM_False);
    for (size_t i = 0; i < n; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);
        if (row->data.function.arg_count != n) return expr_new_symbol(SYM_False);
        for (size_t j = 0; j < n; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol(SYM_False);
            }
        }
    }
    return expr_new_symbol(SYM_True);
}

/* --- DiagonalMatrixQ -----------------------------------------------------
 *
 * `DiagonalMatrixQ[m]` returns True iff every off-(main-diagonal) entry of
 * the matrix `m` is structurally zero; False otherwise.
 *
 * `DiagonalMatrixQ[m, k]` checks the k-th diagonal: positive k selects a
 * superdiagonal above the main diagonal, negative k selects a subdiagonal
 * below it. The matrix is "k-diagonal" iff every entry m[i,j] with
 * j - i != k is zero. Rectangular matrices are supported -- only shape and
 * the entry-zero predicate matter.
 *
 * `DiagonalMatrixQ[m, ..., Tolerance -> t]` relaxes the zero test so that
 * an entry e is considered zero iff `Abs[e] <= t` evaluates to True. This
 * mirrors the behaviour of HermitianMatrixQ / SymmetricMatrixQ.
 *
 * Diagnostics:
 *   - 0 args -> `DiagonalMatrixQ::argt` to stderr, call unevaluated.
 *   - >= 3 positional args (or a non-Rule arg where an option is expected)
 *     -> `DiagonalMatrixQ::nonopt` to stderr, call unevaluated.
 *
 * Shape rejections that return False (rather than unevaluated): non-list
 * input, the empty list `{}`, ragged rows, and any matrix whose entries
 * are themselves Lists (rank >= 3 tensor). `{{}, {}, ...}` (n x 0) is
 * vacuously diagonal and returns True.
 *
 * For symbolic entries the structural zero test is exact: only the
 * literal numeric zeros (Integer 0, Real 0.0, BigInt 0) count as zero.
 * Pure symbols, function calls, and nonzero literals are non-zero, so the
 * predicate is conservative -- it never proves a matrix is diagonal that
 * is not provably so under structural reasoning alone.
 */

static bool diag_entry_is_exact_zero(Expr* e) {
    if (e == NULL) return false;
    if (e->type == EXPR_INTEGER) return e->data.integer == 0;
    if (e->type == EXPR_REAL) return e->data.real == 0.0;
    if (e->type == EXPR_BIGINT) return mpz_sgn(e->data.bigint) == 0;
    return false;
}

static bool diag_entry_under_tolerance(Expr* e, Expr* tol) {
    /* Build LessEqual[Abs[e], tol] and require evaluation to True. */
    Expr** abs_args = malloc(sizeof(Expr*) * 1);
    abs_args[0] = expr_copy(e);
    Expr* abs_call = expr_new_function(expr_new_symbol(SYM_Abs), abs_args, 1);
    free(abs_args);
    Expr* abs_e = eval_and_free(abs_call);

    Expr** le_args = malloc(sizeof(Expr*) * 2);
    le_args[0] = abs_e;
    le_args[1] = expr_copy(tol);
    Expr* le_call = expr_new_function(expr_new_symbol(SYM_LessEqual), le_args, 2);
    free(le_args);
    Expr* le_e = eval_and_free(le_call);
    bool ok = (le_e->type == EXPR_SYMBOL && le_e->data.symbol.name == SYM_True);
    expr_free(le_e);
    return ok;
}

/* Print a Mathematica-compatible argt diagnostic and return NULL so the
 * evaluator leaves the call unevaluated, matching builtin_square_matrix_q
 * / builtin_conjugate's surface behaviour. */
static Expr* diag_emit_argt(size_t argc) {
    fprintf(stderr,
            "DiagonalMatrixQ::argt: DiagonalMatrixQ called with %zu "
            "argument%s; 1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

/* Print a Mathematica-compatible nonopt diagnostic.  We mirror Wolfram's
 * surface text: "Options expected (instead of <expr>) beyond position 2".
 * `bad` is the offending non-Rule expression that broke option parsing. */
static Expr* diag_emit_nonopt(Expr* bad, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "DiagonalMatrixQ::nonopt: Options expected (instead of %s) "
            "beyond position 2 in %s. An option must be a rule or a list "
            "of rules.\n",
            bad_str ? bad_str : "?", call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

Expr* builtin_diagonal_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    ND_MATRIXQ_DELIST(res);
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return diag_emit_argt(0);

    Expr* m = res->data.function.args[0];

    /* Argument layout:
     *   args[0] = matrix m
     *   args[1] = optional integer k (default 0) OR first option Rule.
     *   args[2..] = options.
     */
    int64_t k = 0;
    size_t opt_start = 1;
    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        bool is_rule = (a1->type == EXPR_FUNCTION &&
                        a1->data.function.head->type == EXPR_SYMBOL &&
                        a1->data.function.head->data.symbol.name == SYM_Rule);
        if (is_rule) {
            /* k defaults to 0; options start at position 1. */
            opt_start = 1;
        } else {
            /* Position 1 must be an integer k. */
            if (a1->type != EXPR_INTEGER) {
                /* Non-integer, non-Rule -> treat as a bad option starting
                 * at position 2 (1-indexed). */
                return diag_emit_nonopt(a1, res);
            }
            k = a1->data.integer;
            opt_start = 2;
        }
    }

    /* Parse options. Only Tolerance is recognised; anything else is a
     * nonopt diagnostic (unknown-option rules left over would otherwise
     * silently no-op). */
    Expr* tolerance = NULL;
    Expr* last_bad = NULL;  /* report the LAST non-option encountered */
    for (size_t i = opt_start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        bool is_rule = (opt->type == EXPR_FUNCTION &&
                        opt->data.function.head->type == EXPR_SYMBOL &&
                        opt->data.function.head->data.symbol.name == SYM_Rule &&
                        opt->data.function.arg_count == 2 &&
                        opt->data.function.args[0]->type == EXPR_SYMBOL);
        if (!is_rule) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol.name;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            /* Unknown option name -> treat as a bad option. */
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        return diag_emit_nonopt(last_bad, res);
    }

    /* Validate matrix shape: must be a List of Lists, all rows the same
     * length, no entries that are themselves Lists.  The empty list `{}`
     * (a vector with zero entries) is not a matrix and returns False;
     * `{{}, {}, ...}` (n x 0) is accepted as a vacuous matrix. */
    if (!is_listq(m)) return expr_new_symbol(SYM_False);
    size_t nrows = m->data.function.arg_count;
    if (nrows == 0) return expr_new_symbol(SYM_False);
    size_t ncols = 0;
    for (size_t i = 0; i < nrows; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);
        size_t this_ncols = row->data.function.arg_count;
        if (i == 0) {
            ncols = this_ncols;
        } else if (this_ncols != ncols) {
            return expr_new_symbol(SYM_False);
        }
        for (size_t j = 0; j < this_ncols; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol(SYM_False);
            }
        }
    }

    /* Walk every cell; an entry off the k-th diagonal (j - i != k) must
     * be zero under the chosen predicate.  Entries on the diagonal are
     * unconstrained -- they may be any value, including symbolic. */
    for (size_t i = 0; i < nrows; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = 0; j < ncols; j++) {
            int64_t off = (int64_t)j - (int64_t)i;
            if (off == k) continue;  /* on the k-th diagonal -- ignore */
            Expr* entry = row_i[j];
            bool zero;
            if (tolerance != NULL) {
                zero = diag_entry_under_tolerance(entry, tolerance);
            } else {
                zero = diag_entry_is_exact_zero(entry);
            }
            if (!zero) return expr_new_symbol(SYM_False);
        }
    }

    return expr_new_symbol(SYM_True);
}

/* --- UpperTriangularMatrixQ --------------------------------------------
 *
 * `UpperTriangularMatrixQ[m]` returns True iff every entry of `m` strictly
 * below the main diagonal is zero; False otherwise.
 *
 * `UpperTriangularMatrixQ[m, k]` shifts the cut-off to the k-th diagonal:
 * every entry m[i,j] with `j - i < k` must be zero.  Equivalently, all
 * nonzero entries are confined to the region on or above the k-th
 * diagonal.  Positive k selects a superdiagonal above the main diagonal
 * (so the test becomes stricter); negative k selects a subdiagonal below
 * it (the test becomes more permissive).  Rectangular matrices are
 * supported.
 *
 * `UpperTriangularMatrixQ[m, ..., Tolerance -> t]` relaxes the zero test
 * so that an entry e is considered zero iff `Abs[e] <= t` evaluates to
 * True.  Mirrors the surface API of DiagonalMatrixQ / HermitianMatrixQ /
 * SymmetricMatrixQ.
 *
 * Diagnostics (Mathematica-compatible, to stderr):
 *   - 0 args -> `UpperTriangularMatrixQ::argt`, call left unevaluated.
 *   - >= 3 positional args (or a non-Rule arg where an option is
 *     expected) -> `UpperTriangularMatrixQ::nonopt`, call left
 *     unevaluated.
 *
 * Shape rejections that return False (rather than unevaluated): non-list
 * input, the empty list `{}`, ragged rows, and any matrix whose entries
 * are themselves Lists (rank >= 3 tensor).  `{{}, {}, ...}` (n x 0) is
 * vacuously upper-triangular and returns True.
 *
 * For symbolic entries the structural zero test is exact: only literal
 * numeric zeros (Integer 0, Real 0.0, BigInt 0) count as zero.  Symbolic
 * sub-diagonal entries cause the matrix to be rejected, so the predicate
 * is conservative.
 */

static Expr* utri_emit_argt(size_t argc) {
    fprintf(stderr,
            "UpperTriangularMatrixQ::argt: UpperTriangularMatrixQ called "
            "with %zu argument%s; 1 or 2 arguments are expected.\n",
            argc, argc == 1 ? "" : "s");
    return NULL;
}

static Expr* utri_emit_nonopt(Expr* bad, Expr* res) {
    char* bad_str = expr_to_string(bad);
    char* call_str = expr_to_string(res);
    fprintf(stderr,
            "UpperTriangularMatrixQ::nonopt: Options expected (instead of "
            "%s) beyond position 2 in %s. An option must be a rule or a "
            "list of rules.\n",
            bad_str ? bad_str : "?", call_str ? call_str : "?");
    free(bad_str);
    free(call_str);
    return NULL;
}

Expr* builtin_upper_triangular_matrix_q(Expr* res) {
    if (res->type != EXPR_FUNCTION) return NULL;
    ND_MATRIXQ_DELIST(res);
    size_t argc = res->data.function.arg_count;
    if (argc == 0) return utri_emit_argt(0);

    Expr* m = res->data.function.args[0];

    /* Argument layout, mirroring DiagonalMatrixQ:
     *   args[0]    = matrix m
     *   args[1]    = optional integer k (default 0) OR first option Rule.
     *   args[2..]  = options.
     */
    int64_t k = 0;
    size_t opt_start = 1;
    if (argc >= 2) {
        Expr* a1 = res->data.function.args[1];
        bool is_rule = (a1->type == EXPR_FUNCTION &&
                        a1->data.function.head->type == EXPR_SYMBOL &&
                        a1->data.function.head->data.symbol.name == SYM_Rule);
        if (is_rule) {
            opt_start = 1;
        } else {
            if (a1->type != EXPR_INTEGER) {
                return utri_emit_nonopt(a1, res);
            }
            k = a1->data.integer;
            opt_start = 2;
        }
    }

    /* Parse options.  Only Tolerance is recognised. */
    Expr* tolerance = NULL;
    Expr* last_bad = NULL;
    for (size_t i = opt_start; i < argc; i++) {
        Expr* opt = res->data.function.args[i];
        bool is_rule = (opt->type == EXPR_FUNCTION &&
                        opt->data.function.head->type == EXPR_SYMBOL &&
                        opt->data.function.head->data.symbol.name == SYM_Rule &&
                        opt->data.function.arg_count == 2 &&
                        opt->data.function.args[0]->type == EXPR_SYMBOL);
        if (!is_rule) {
            last_bad = opt;
            continue;
        }
        const char* name = opt->data.function.args[0]->data.symbol.name;
        Expr* val = opt->data.function.args[1];
        if (name == SYM_Tolerance) {
            if (!(val->type == EXPR_SYMBOL && val->data.symbol.name == SYM_Automatic)) {
                tolerance = val;
            }
        } else {
            last_bad = opt;
        }
    }
    if (last_bad != NULL) {
        return utri_emit_nonopt(last_bad, res);
    }

    /* Validate matrix shape -- same gate as DiagonalMatrixQ. */
    if (!is_listq(m)) return expr_new_symbol(SYM_False);
    size_t nrows = m->data.function.arg_count;
    if (nrows == 0) return expr_new_symbol(SYM_False);
    size_t ncols = 0;
    for (size_t i = 0; i < nrows; i++) {
        Expr* row = m->data.function.args[i];
        if (!is_listq(row)) return expr_new_symbol(SYM_False);
        size_t this_ncols = row->data.function.arg_count;
        if (i == 0) {
            ncols = this_ncols;
        } else if (this_ncols != ncols) {
            return expr_new_symbol(SYM_False);
        }
        for (size_t j = 0; j < this_ncols; j++) {
            if (is_listq(row->data.function.args[j])) {
                return expr_new_symbol(SYM_False);
            }
        }
    }

    /* Every entry strictly below the k-th diagonal (j - i < k) must be
     * zero under the chosen predicate.  Entries on or above the k-th
     * diagonal are unconstrained. */
    for (size_t i = 0; i < nrows; i++) {
        Expr** row_i = m->data.function.args[i]->data.function.args;
        for (size_t j = 0; j < ncols; j++) {
            int64_t off = (int64_t)j - (int64_t)i;
            if (off >= k) continue;  /* on or above k-th diagonal -- ignore */
            Expr* entry = row_i[j];
            bool zero;
            if (tolerance != NULL) {
                zero = diag_entry_under_tolerance(entry, tolerance);
            } else {
                zero = diag_entry_is_exact_zero(entry);
            }
            if (!zero) return expr_new_symbol(SYM_False);
        }
    }

    return expr_new_symbol(SYM_True);
}
